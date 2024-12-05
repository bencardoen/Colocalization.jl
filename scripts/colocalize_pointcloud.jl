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

f1 = "/home/bcardoen/cedar_data/mattest/two_channel_mat/abbelight/CavPTRF_2_2_Cav_647_merged_threshold_0_alpha_3.mat"
f2 = "//home/bcardoen/cedar_data/mattest/two_channel_mat/abbelight/CavPTRF_2_2_PTRF_680_merged_threshold_0_alpha_3.mat"

# dfx = coloc_srn(f1, f2)

### Array cluster ID to class ID
# cluster_to_class = data["clstClass"]

### 3D points 
# points=data["DatFiltered"]

# clustercount = size(cluster_to_class, 2)

# p2c = data["point2cluster"]

# c2p = data["clustMembsCell"]

# c28 = Int.(c2p[i])
# p28 = points[c28, :][1,:,:]

# ctr = mean(p28, dims=1)


# p2c1, c2p1, clustercount1, points1, cluster_to_class1 = load_SRN(f1)
# ctrs1, segs1, L1 = compute_centroids(points1, c2p1)

# channel_ctr1 = mean(ctrs1, dims=1)[:]


# p2c2, c2p2, clustercount2, points2, cluster_to_class2 = load_SRN(f2)
# ctrs2, segs2, L2 = compute_centroids(points2, c2p2)

# channel_ctr2 = mean(ctrs2, dims=1)[:]

# df12 = report_distances(ctrs1, ctrs2, channel_ctr1)
# df21 = report_distances(ctrs2, ctrs1, channel_ctr2)
# df12[!,:channel].=1
# df21[!,:channel].=2
# dfx=vcat([df12, df21]...)


function report_distances(ctrs1, ctrs2, channel_ctr1)
    df = DataFrame(clusterid1 = Int[], centroid_x=Float64[], centroid_y=Float64[], centroid_z=Float64[], 
    distance_1 = Float64[], distance_2 = Float64[], distance_3 = Float64[], distance_4 = Float64[], distance_5 = Float64[], 
    nearest_1 = Int[], nearest_2 = Int[], nearest_3 = Int[], nearest_4 = Int[], nearest_5 = Int[],
    channel_centroid_x=Float64[], channel_centroid_y=Float64[], 
    channel_centroid_z=Float64[], distance_to_centroid=Float64[])
    for (i, c1) in enumerate(eachrow(ctrs1))
        d1 = (ctrs2 .- reshape(c1, (1,3))).^2
        d1s = sqrt.(sum(d1, dims=2))[:]
        sd = sortperm(d1s)
        near = sd[1:5]
        neardistance = d1s[sd][1:5]
        push!(df, [i, c1[1], c1[2], c1[3], neardistance..., near..., channel_ctr1[1], channel_ctr1[2], channel_ctr1[3], sqrt(sum((c1 .- channel_ctr1).^2))])
    end
    return df
end

function compute_centroids(ps, cs, minsize=5)
    CN = size(cs, 1)
    centroids = zeros(CN, 3)
    segments = Dict()
    lengths = zeros(CN)
    for i in 1:CN
        ip = Int.(cs[i])
        ni = length(ip)
        lengths[i] = ni
        if ni < minsize
            @warn "Have < $minsize points for cluster $i"
        end
        if ni < 2
            ip = ip
        else
            ip = ip[:]
        end
        ipcoords = ps[ip, :]
        segments[i] = ipcoords
        centroids[i, :] .= mean(ipcoords, dims=1)[1,:]
    end
    return centroids, segments, lengths
end

function load_SRN(matfile)
    data = matread(matfile)
    cluster_to_class = data["clstClass"]
    points=data["DatFiltered"]
    clustercount = size(cluster_to_class, 2)
    p2c = data["point2cluster"]
    c2p = data["clustMembsCell"]
    return p2c, c2p, clustercount, points, cluster_to_class
end

function coloc_srn(f1, f2)
    p2c1, c2p1, clustercount1, points1, cluster_to_class1 = load_SRN(f1)
    ctrs1, segs1, L1 = compute_centroids(points1, c2p1)
    channel_ctr1 = mean(ctrs1, dims=1)[:]
    p2c2, c2p2, clustercount2, points2, cluster_to_class2 = load_SRN(f2)
    ctrs2, segs2, L2 = compute_centroids(points2, c2p2)
    channel_ctr2 = mean(ctrs2, dims=1)[:]
    df12 = report_distances(ctrs1, ctrs2, channel_ctr1)
    df21 = report_distances(ctrs2, ctrs1, channel_ctr2)
    df12[!,:channel].=1
    df12[!,:channelfile].=f1
    df21[!,:channel].=2
    df21[!,:channelfile].=f2
    dfx=vcat([df12, df21]...)
    return dfx
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--first", "-f"
            help = "Filename of first channel"
            arg_type = String
            required = true
		"--second", "-s"
            help = "Filename of second channel"
            arg_type = String
            required = true
		"--outdir", "-o"
            help = "output folder"
            arg_type = String
            required = true
    end

    return parse_args(s)
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
	first = parsed_args["first"]
	second = parsed_args["second"]
	outdir = parsed_args["outdir"]
	if ! all(isfile.([first, second]))
		@error "Missing files $first $second"
		return
	end
	df = coloc_srn(first, second)
	CSV.write("results.csv", dfx)
	@info "Done"
end


runcoloc()