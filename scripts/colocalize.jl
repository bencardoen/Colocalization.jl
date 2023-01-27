
# Copyright (C) 2018-2023 Ben Cardoen bcardoen@sfu.ca
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
using CSV, Statistics
import Glob
using Logging
using Dates
using ArgParse
using LoggingExtras
using SmlmTools
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
            help = "output folder"
            arg_type = String
            required = true
		"--segment", "-g"
            help = "If active, segment the 2D images before colocalization. Not all coloc methods work without segmentation."
            action = :store_true
            default = false
		"--windowsize", "-w"
			help = "Windowsize used in colocalization (default = 3, should be odd >= 3)"
			arg_type = Int64
			default = 3
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
	C1 = Images.load(first)
	C2 = Images.load(second)
	if parsed_args["segment"]
		@info "Segmenting"
		C1 = segment(C1)
		C2 = segment(C2)
	end
	results = colocalize_all(C1, C2; windowsize=parsed_args["windowsize"])
	for k in keys(results)
		Images.save(joinpath(outdir, "$k.tif"), N0f16.(nmz(abs.(results[k]))))
	end
	df=summarize_colocalization(results, first, second)
	CSV.write(joinpath(outdir, "colocalization_results.csv"), df)
	@info "Done"
end


runalign()
