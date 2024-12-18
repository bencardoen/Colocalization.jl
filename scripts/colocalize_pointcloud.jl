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
            required = true
		"--second", "-s"
            help = "Filename of second channel"
            arg_type = String
            required = true
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
	df = coloc_srn(first, second; SRN_minpoints=parsed_args["mincluster"])
    @info "Done... saving to CSV: $(joinpath(outdir, "$(outname).csv"))"
	CSV.write(joinpath(outdir, "$(outname).csv"), df)
	@info "Done"
end


runcoloc()