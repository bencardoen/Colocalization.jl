# Copyright (C) 2023 Ben Cardoen bcardoen@sfu.ca
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


module Colocalization

using Statistics
using DataFrames
using ProgressMeter
using Images
using Distributions
using ImageFiltering

export describe_array, colocalize_all, colocalize, segment, tomask, aszero, summarize_colocalization


"""
	describe_array(xs)

	Returns μ, md, σ, min, max, q1, q3, iqr, q95, q99, nans for `xs`.

"""
function describe_array(xs)
	xs = Float64.(xs)
    X, Y = size(xs)
    nans = length(xs[isnan.(xs)])
    if nans == X*Y
        return NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, nans
    end
    xps = xs[.! isnan.(xs)]
	m, M = minimum(xps), maximum(xps)
	md = quantile(xps, [.5])[1]
	μ = mean(xps)
	σ = std(xps)
    q1, q3, q95, q99 = quantile(xps, [0.25, 0.75, 0.95, 0.99])
    iqr = (q3-q1)/2
    return μ, md, σ, m, M, q1, q3, iqr, q95, q99, nans
end

"""
    tomask(img)

	Utility function to binarize argument
"""
function tomask(img)
    c = copy(img)
    c[c .> 0] .= 1
    return Images.Gray{Images.N0f8}.(c)
end

"""
	aszero(xs)

	zeros(eltype(xs), sizes(xs))
"""
function aszero(xs)
    return zeros(eltype(xs), size(xs))
end

"""
	colocalize_all(xs, ys; windowsize=3)

	Apply all coloc metrics. Returns a dictionary, so results are stored in results[metric].
	Runs multithreaded.
	```julia
	xs = rand(10, 10)
	ys = rand(10, 10)
	res =colocalize_all(xs, ys)
	correlation_values = res["pearson"]
	μ, md, σ, min, max, q1, q3, iqr, q95, q99, nans = describe_map(correlation_values)
	```
"""
function colocalize_all(xs, ys; windowsize=3)
    mt = keys(Colocalization.metrics) |> collect
    res = Dict([(k, zeros(size(xs))) for k in mt])
    for metric in mt
        res[metric] = colocalize(xs, ys; metric=metric, windowsize=windowsize)
    end
    return res
end

function segment(img)
    rimg = copy(img)
    t=otsu_threshold(rimg)
    rimg[rimg .< t] .=0
    return rimg
end

"""
	colocalize(xs, ys; metric="pearson", windowsize=3)

	Compute a colocalization metric `metric` on a window size of `windowsize^2`.

	See `metrics` for metrics to use.
"""
function colocalize(_xs, _ys; metric="pearson", windowsize=3)
	xs = Float64.(_xs)
	ys = Float64.(_ys)
	@info "Coloc with window $windowsize  and metric $metric for input $(size(xs))"
    X, Y = size(xs)
	if size(xs)!=size(ys)
		@error "Dimensions $(size(xs)) $(size(ys))"
		throw(ArgumentError("Dimensions $(size(xs)) $(size(ys))"))
	end
	if !((windowsize %2 == 1) || (windowsize >= 3))
		@error "Invalid windowsize $windowsize , should be >=3 and odd"
		throw(ArgumentError("Invalid windowsize $windowsize , should be >=3 and odd"))
	end
    k = Int((windowsize-1)/2)
	if !(metric ∈ keys(metrics))
		@error "Invalid metric: supported metrics are $(keys(metrics))"
		throw(ArgumentError("Invalid metric $metric"))
	end
    # metrics = Dict([("pearson", coloc_pearson), ("spearman", coloc_spearman), ("jaccard", coloc_jaccard), ("manders", coloc_manders), ("m1", coloc_m1), ("sorensen", coloc_sorensen), ("m2", coloc_m2)])
    mf = metrics[metric]
    result = aszero(xs)
	N = (X-k*2)
    @showprogress for yi in k+1:(Y-k)
			for xi in k+1:(X-k)
    			result[xi, yi] = nanz(mf(xs[xi-k:xi+k, yi-k:yi+k], ys[xi-k:xi+k, yi-k:yi+k]))
			end
    end
    return result
end

function summarize_colocalization(results, f1, f2)
	cols = ["μ", "md", "σ", "min", "max", "q1", "q3", "iqr", "q95", "q99", "nans"]
    vals = zeros(length(keys(Colocalization.metrics)), length(cols))
    meta = []
    for (i,m) in enumerate(keys(Colocalization.metrics))
        vals[i,:] .= describe_array(results[m])
        push!(meta, [m])
    end
    df=DataFrame(vals, cols)
    df[!, :channel_1] .= f1
    df[!, :channel_2] .= f2
    df[!, :metric] .= meta
    return df
end

function nanz(x)
	isnan(x) ? 0 : x
end

function coloc_pearson(xs, ys)
        cor(xs[:], ys[:])
end

function coloc_spearman(xs, ys)
    return spear(xs[:], ys[:])[1]
end

function coloc_manders(xs, ys)
    xs = xs[:]
    ys = ys[:]
    return sum(xs .* ys) / sqrt( sum(xs.*xs) + sum(ys.*ys) )
end

function coloc_jaccard(xs, ys)
    xs = xs[:]
    ys = ys[:]
    sum((xs .> 0) .& (ys .> 0)) / sum((xs .> 0) .| (ys .> 0))
end

function coloc_sorensen(xs, ys)
    xs = xs[:]
    ys = ys[:]
    j = coloc_jaccard(xs, ys)
    return (2*j)/(1+j)
end

function coloc_m1(xs, ys)
    xs = xs[:]
    ys = ys[:]
    sum(xs[ys .> 0])/sum(xs)
end

function coloc_m2(xs, ys)
    xs = xs[:]
    ys = ys[:]
    sum(ys[xs .> 0])/sum(ys)
end

metrics = Dict([("pearson", coloc_pearson), ("spearman", coloc_spearman), ("jaccard", coloc_jaccard), ("manders", coloc_manders), ("m1", coloc_m1), ("sorensen", coloc_sorensen), ("m2", coloc_m2)])



"""
	spear(left, right)
	Return the spearman correlation between arguments, and the significance (z-test, t-test)

	Reuse with permission from https://github.com/bencardoen/SubPrecisionContactDetection.jl
"""
function spear(view1, view2)
    N = length(view1)
    srt1 = sortperm(view1[:])
    srt2 = sortperm(view2[:])
    r = cor(srt1, srt2)
	z, t = computezcorr(r, N)
    return r, z, t
end

function computezcorr(r, N)
	# Fisher transform
	fr = atanh(r)
    z = sqrt((N-3)/1.06) * fr
    t = r * (sqrt((N-2)/(1-r^2)))
	return z, t
end



end
