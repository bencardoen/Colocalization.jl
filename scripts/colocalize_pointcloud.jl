# Copyright (C) 2024 Ben Cardoen bcardoen@sfu.ca
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
using Pkg; Pkg.activate(".")
using MAT, DataFrames, CSV
using ArgParse
using CSV, Statistics
import Glob
using Logging
using Dates
using ArgParse
using LoggingExtras
using Colocalization
using Images


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--first", "-f"
            help = "Filename of first channel"
            arg_type = String
            required = false
		"--second", "-s"
            help = "Filename of second channel"
            arg_type = String
            required = false
        "--interaction", "-i"
            help = "Interaction factor, default 1. <, select for tighter interaction, > select for long range interaction"
            arg_type = Float64
            required = false
            default = 1.0
        "--patternmatch", "-p"
            help = "Search for matching pairs of files in this folder"
            arg_type = String
            required = false
            default = ""
        "--pattern", "-r"
            help = "Pattern to use to find pairs of files"
            arg_type = String
            required = false
            default = "CavPTRF_1_1_[A-Z, a-z]+_[0-9]+_merged_threshold_[0-9]+_alpha_[0-9]+.mat"    
        "--outdir", "-o"
            help = "output folder, default to current directory"
            arg_type = String
            default = ""
            required = false
        "--mincluster", "-m"
            help = "Minimum size of cluster (nr of points). Default to 4, NEVER set < 4. Should correspond to SRN choice."
            arg_type = Int
            default = 4
    end

    return parse_args(s)
end


function find_pairs(dir, pattern)
    files = readdir(dir)
    filtered = [f for f in files if !isnothing(match(pattern, f))]
    pairs = Dict()
    for f in filtered
        try 
            parts = split(split(f, ".")[1], '_')
            t, a = parts[end-2], parts[end]
            k = (t, a)
            if k in keys(pairs)
                push!(pairs[k], joinpath(dir, f))
            else
                pairs[k] = [joinpath(dir, f)]
            end
        catch  e 
            @warn "Parsing $f failed due to $e"
        end 
    end
    @info "Checking pairs"
    for k in keys(pairs)
        if length(pairs[k]) != 2
            @warn "Not a pair for $(k) $(pairs[k])"
            delete!(pairs, k)
        end
    end
    @info "Found a total of $(length(keys(pairs)))"
    return pairs
end


function pairwise(parsed_args)
    first = parsed_args["first"]
	second = parsed_args["second"]
	outdir = parsed_args["outdir"]
    F= parsed_args["interaction"]
    if F <= 0
        @error "Invalid interaction value $F"
        F = 1.0
    end
    if outdir == ""
        @info "No output directory given, working in current directory $(pwd())"
        outdir = pwd()
    end
	if ! all(isfile.([first, second]))
		@error "Missing files $first $second"
		return
	end
    f = splitext(splitpath(first)[end])[1]
    s = splitext(splitpath(second)[end])[1]
    outname = "$(f)___$(s)"
    if length(outname) > 255
        @warn "Output name probably too long ..."
    end
    @info "Starting Colocalization"
	dfx, data = coloc_srn(first, second; SRN_minpoints=parsed_args["mincluster"])
    @info "Done... saving to CSV: $(joinpath(outdir, "$(outname).csv"))"
	CSV.write(joinpath(outdir, "$(outname).csv"), dfx)
    @info "Writing pointclouds"
    segs1=data["segs1"]
    c1=dfx[dfx.channel .== 1, :]
    # write_clusters(segs1, "channel_1", outdir)
    int_segs = filter_interacting(c1, segs1, Inf)
    write_clusters(int_segs, "channel_1", outdir)
    int_segs = filter_interacting(c1, segs1, F)
    write_clusters(int_segs, "channel_1_interacting", outdir)
    pts=vcat([int_segs[k] for k in keys(int_segs)]...)
    points_to_vtu(pts, joinpath(outdir, "c1_interacting"))
    points_to_vtu(data["points1"], joinpath(outdir, "c1_all"))

    segs2=data["segs2"]
    c2=dfx[dfx.channel .== 2, :]
    # write_clusters(segs2, "channel_2", outdir)
    int_segs = filter_interacting(c2, segs2, Inf)
    write_clusters(int_segs, "channel_2", outdir)
    int_segs = filter_interacting(c2, segs2, F)
    write_clusters(int_segs, "channel_2_interacting", outdir)
    pts=vcat([int_segs[k] for k in keys(int_segs)]...)
    points_to_vtu(pts, joinpath(outdir, "c2_interacting"))
    points_to_vtu(data["points2"], joinpath(outdir, "c2_all"))
	@info "Done"
end

function runcoloc()
    date_format = "yyyy-mm-dd HH:MM:SS"
    timestamp_logger(logger) = TransformerLogger(logger) do log
      merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
    end
    ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        @info "  $arg  =>  $val"
    end
	@info "Finished"
    ### If pattern is given, the find pairs and proceed 
    ### If not, then as usual
    pmatch = parsed_args["patternmatch"]
    if pmatch == ""
        pairwise(parsed_args)
    else
        pattern = parsed_args["pattern"]
        if !isdir(pmatch)
            @error "$pmatch is not a directory, quitting."
            exit(-1)
        end
        pairs = find_pairs(pmatch, Regex(pattern))
        args = copy(parsed_args)
        outdir = parsed_args["outdir"]
        if outdir == ""
            @info "Outdir is empty, setting to current directory."
            outdir = pwd()
        end
        for k in keys(pairs)
            files = sort(pairs[k])
            #  Set f, s 
            # Set outdir
            t, a = k
            args = copy(parsed_args)
            @info "$t $a $files"
            args["first"] = files[1]
            args["second"] = files[2]
            args["outdir"] = joinpath(outdir, "$(t)_$(a)")
            mkpath(args["outdir"])
            @info "Saving to $(args["outdir"])"
            pairwise(args)
        end
    end
	# first = parsed_args["first"]
	# second = parsed_args["second"]
	# outdir = parsed_args["outdir"]
    # F= parsed_args["interaction"]
    # if F <= 0
    #     @error "Invalid interaction value $F"
    #     F = 1.0
    # end
    # if outdir == ""
    #     @info "No output directory given, working in current directory $(pwd())"
    #     outdir = pwd()
    # end
	# if ! all(isfile.([first, second]))
	# 	@error "Missing files $first $second"
	# 	return
	# end
    # f = splitext(splitpath(first)[end])[1]
    # s = splitext(splitpath(second)[end])[1]
    # outname = "$(f)___$(s)"
    # if length(outname) > 255
    #     @warn "Output name probably too long ..."
    # end
    # @info "Starting Colocalization"
	# dfx, data = coloc_srn(first, second; SRN_minpoints=parsed_args["mincluster"])
    # @info "Done... saving to CSV: $(joinpath(outdir, "$(outname).csv"))"
	# CSV.write(joinpath(outdir, "$(outname).csv"), dfx)
    # @info "Writing pointclouds"
    # segs1=data["segs1"]
    # c1=dfx[dfx.channel .== 1, :]
    # # write_clusters(segs1, "channel_1", outdir)
    # int_segs = filter_interacting(c1, segs1, Inf)
    # write_clusters(int_segs, "channel_1", outdir)
    # int_segs = filter_interacting(c1, segs1, F)
    # write_clusters(int_segs, "channel_1_interacting", outdir)
    # pts=vcat([int_segs[k] for k in keys(int_segs)]...)
    # points_to_vtu(pts, joinpath(outdir, "c1_interacting"))
    # points_to_vtu(data["points1"], joinpath(outdir, "c1_all"))

    # segs2=data["segs2"]
    # c2=dfx[dfx.channel .== 2, :]
    # # write_clusters(segs2, "channel_2", outdir)
    # int_segs = filter_interacting(c2, segs2, Inf)
    # write_clusters(int_segs, "channel_2", outdir)
    # int_segs = filter_interacting(c2, segs2, F)
    # write_clusters(int_segs, "channel_2_interacting", outdir)
    # pts=vcat([int_segs[k] for k in keys(int_segs)]...)
    # points_to_vtu(pts, joinpath(outdir, "c2_interacting"))
    # points_to_vtu(data["points2"], joinpath(outdir, "c2_all"))
	# @info "Done"
end


runcoloc()