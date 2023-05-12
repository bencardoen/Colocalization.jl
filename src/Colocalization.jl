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
using ImageMorphology

export describe_array, intersection_mask, haussdorff_max, haussdorff_mean, haussdorff_distance, filter_projection, object_stats, union_mask, union_distance_mask, colocalize_nowindow, colocalize_all, colocalize, segment, tomask, aszero, summarize_colocalization, list_metrics, metrics_iterator


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
	filter_projection(image, w=1, z=0)

	Compute a 2-stage filtered mask based on image. w is the window size, z is the number of standard deviations to use as a threshold.
	This removes isolated localizations of low intensity.
"""
function filter_projection(im,w=1, z=0)
    # Source https://github.com/bencardoen/SPECHT.jl/blob/main/src/SPECHT.jl
    μ , _,σ = describe_array(im)[1:3]
    i1 = copy(im)
    i1[i1 .< μ + z*σ] .= 0
    mg= Images.mapwindow(Statistics.median, i1, [w*2+1 for _ in 1:length(size(im))])
    mg[isnan.(mg)] .= zero(eltype(im))
    mg= Images.mapwindow(Statistics.median, mg, [w*2+1 for _ in 1:length(size(im))])
    mg[isnan.(mg)] .= zero(eltype(im))
    return tomask(im) .* tomask(mg) # Make sure we don't include FP
end

"""
	haussdorff_distance(xs, ys, aggregate=maximum)

	Returns, for two segmented masks, the Haussdorff distance (per object).

	Agg=maximum is the default, returns the maximum symmetric Euclidean distance of X to Y
"""
function _haussdorff_distance(A, B; agg=maximum)
    distA = distance_transform(feature_transform(A .> 0));
    # @info eltype(distA)
    d = zeros(size(B))
    d .== Inf
    if sum(A) == 0
        return d
    end
    if sum(B) == 0
        return d
    end
    idsB = component_indices(label_components(B))[2:end]
    mBA = zeros(Float64, size(B))
    for ind in idsB
        mBA[ind] .= agg(distA[ind])
    end
    distB = distance_transform(feature_transform(B .> 0));
    idsA = component_indices(label_components(A))[2:end]
    mAB = zeros(size(B))
    for ind in idsA
        mAB[ind] .= agg(distB[ind])
    end 
    return mBA, mAB, max.(mBA, mAB)
end


"""
	haussdorff_distance(xs, ys, aggregate=maximum)

	Returns, for two segmented masks, the Haussdorff distance (per object).

	Agg=maximum is the default, returns the maximum symmetric Euclidean distance of X to Y
"""
function haussdorff_distance(A, B; agg=maximum)
	#@info agg
    if sum(A) == 0
        # @warn "No pixels in A"
        return 0
    end
    if sum(B) == 0
        # @warn "No pixels in B"
        return 0
    end
    distA = distance_transform(feature_transform(A .> 0));
    dAB = agg(distA[B .> 0])
    distB = distance_transform(feature_transform(B .> 0));
    dBA = agg(distB[A .> 0]) 
    return max(dAB, dBA)
end

function _rlap(lap)
	# Reused from SPCHT
    return _rlap!(copy(lap))
end

function _rlap!(lap)
	# Reused from SPECHT
	z = zero(eltype(lap))
	_f = x ->  x < z ? -x : z
	return map!(_f, lap, lap)
end

function _find_th(imgl, Z, auto=false, scale=1, edgemask=nothing)
	negl = _rlap(imgl)
	if ! isnothing(edgemask)
		@debug "Using edgemask"
		@assert size(edgemask) == size(imgl)
		negl[edgemask .> 0 ] .= 0
	end
    Zp = Z
	avals = Float64.(negl[:])
    if auto == true
		kurt = Distributions.kurtosis(avals)
        Zp = kurt^(.25)/scale
        @debug "Using self scaling Z of $(Zp)"
    end
    gem, ges = _gmsm(avals)
    Tg = gem * ges^Zp
    return Tg, negl
end


function _gmsm(vs)
    # Reused from ERGO with permission
	vals = copy(vs)
    vals = vals[vals .> zero(eltype(vals))]
    N = length(vals)
    li = [log(i) for i in vals]
    gm = exp(sum(li) / N)
    gs = exp(sqrt(sum([log(i/gm)^2 for i in vals])/N))
    return gm, gs
end


"""
	haussdorff_max(xs, ys)

	The reference implementation for the Haussdorff distance, returns the maximum symmetric Euclidean distance of X to Y
	See https://en.wikipedia.org/wiki/Hausdorff_distance
"""
function haussdorff_max(xs, ys)
	# @info "max"
	return haussdorff_distance(xs, ys, agg=maximum)
end


"""
	haussdorff_mean(xs, ys)

	A variant of the Haussdorff distance, returns the mean symmetric Euclidean distance of X to Y
	See https://en.wikipedia.org/wiki/Hausdorff_distance
"""
function haussdorff_mean(xs, ys)
	# @info "mean"
	return haussdorff_distance(xs, ys, agg=mean)
end

"""
	union_distance_mask(xs, ys, distance)

	Finds the connected components in X that are <= distance to any component in Y and vice versa

	Returns mask_x, mask_y
"""
function union_distance_mask(xs, ys, distance)
	xm=Int8.(Colocalization.tomask(xs))
	ym=Int8.(Colocalization.tomask(ys))
	distance_map_x = Images.distance_transform(Images.feature_transform(Bool.(xm)))
	distance_map_y = Images.distance_transform(Images.feature_transform(Bool.(ym)))
	rx = Colocalization.aszero(xs)
	ry = Colocalization.aszero(ys)
	ccx = label_components(Colocalization.tomask(xm))
	ccy = label_components(Colocalization.tomask(ym))
	for ind in component_indices(ccx)[2:end]
		if minimum(distance_map_y[ind]) <= distance
			rx[ind] .= 1
		end
	end
	for ind in component_indices(ccy)[2:end]
		if minimum(distance_map_x[ind]) <= distance
			ry[ind] .= 1
		end
	end
	return rx, ry
end

"""
	intersection_mask(xs, ys)

	mask where x .> 0 and y .> 0
"""
function intersection_mask(xs, ys)
	xm=Colocalization.tomask(xs)
	ym=Colocalization.tomask(ys)
	return xm .* ym
end


"""
	union_mask(xs, ys)

	mask where any component that has a non-zero intersection is set to 1.
"""
function union_mask(xs, ys)
	xm=Int8.(Colocalization.tomask(xs))
	ym=Int8.(Colocalization.tomask(ys))
	joint = xm .+ ym
	res = Colocalization.aszero(xs)
	ccs = label_components(Colocalization.tomask(joint))
	for ind in component_indices(ccs)
		if maximum(joint[ind]) == 2
			res[ind] .= 1
		end
	end
	return res
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
function colocalize_all(xs, ys; windowsize=3, metric=nothing)
    mt = keys(Colocalization.metrics) |> collect
	if !isnothing(metric)
		@info "Filtering on $metric"
		filter!(x -> x == metric, mt)
	end
    res = Dict([(k, zeros(size(xs))) for k in mt])
    for met in mt
        res[met] = colocalize(xs, ys; metric=met, windowsize=windowsize)
    end
    return res
end

"""
	segment(image, scale=1.0)

	Perform basic otsu thresholding. Scale argument allows you to include more (>1) background or less (<1)
"""
function segment(img, scale=1.0; method="otsu")
	if method == "otsu"
		rimg = copy(img)
		t=otsu_threshold(rimg) * scale
		rimg[rimg .< t] .=0
		rimg[rimg .> 0] .=1
		return rimg
	end
	if method == "specht"
		return segment_specht(img, scale)
	end
	raise(ArgumentError("Invalid method $method"))
end

function _glg(img, sigmas)
    @assert(length(sigmas)==2)
    imgg = img
    if sigmas[1] != 0
        imgg = ImageFiltering.imfilter(img, ImageFiltering.Kernel.gaussian(sigmas[1]));
    end
    # imgg = img
    imgl = ImageFiltering.imfilter(imgg, ImageFiltering.Kernel.Laplacian());
    imglg = imgl
    if sigmas[2] != 0
        imglg = ImageFiltering.imfilter(imgl, ImageFiltering.Kernel.gaussian(sigmas[2]));
    end
    return img, imgg, imgl, imglg
end


function segment_specht(img, prc=2, sigmas=[2,2])
    #Reused with permission from https://github.com/bencardoen/SPECHT.jl/blob/main/src/SPECHT.jl
    _, _, _, imgl = _glg(img, sigmas); 
    Tg, neg_glog = _find_th(imgl, 0, true, prc); 
    ngl = _rlap(imgl)
    i2 = copy(img)
    i2[ngl .< Tg] .= zero(eltype(img))
    i2[ngl .>= Tg] .= oneunit(eltype(img))
    i2[1:3,:].=0
    i2[:,1:3].=0
    i2[end-2:end,:].=0
    i2[:,end-2:end].=0
    return i2
end

function object_stats(i1, i2, context)
    m1= tomask(i1)
    m2= tomask(i2)
    c1 = Images.label_components(m1)
    c2 = Images.label_components(m2)
    disc1_c2, o12,  d1map = compute_distances_cc_to_mask(c1, m2)
    disc2_c1, o21, d2map  = compute_distances_cc_to_mask(c2, m1)
    c1stats = describe_cc(c1, i1)
    c2stats = describe_cc(c2, i2)
    df1 = DataFrame(channel=1, overlap=o12, distance_to_nearest=disc1_c2, area=Images.component_lengths(c1)[2:end], mean=c1stats[:,1], std=c1stats[:,2])
    df2 = DataFrame(channel=2, overlap=o21, distance_to_nearest=disc2_c1, area=Images.component_lengths(c2)[2:end], mean=c2stats[:,1], std=c2stats[:,2])
    return vcat(df1, df2), d1map, d2map
end


function describe_cc(ccs, img)
    N = maximum(ccs)
    dis = zeros(N, 2)
    for (i, ci) in enumerate(component_indices(ccs)[2:end])
        im = img[ci]
        dis[i, :] .= mean(im), std(im)
    end
    return dis
end



function compute_distances_cc_to_mask(from_cc, to_mask)
    # reused from SPECHT https://github.com/bencardoen/SPECHT.jl/blob/main/src/SPECHT.jl
	@assert(size(from_cc) == size(to_mask))
    dismap = Images.distance_transform(Images.feature_transform(Bool.(to_mask)))
    N = Int(maximum(from_cc))
    @info N
    dis = zeros(N)
	overlap = zeros(N)
	distmap = zeros(Float64, size(to_mask))
    ind = Images.component_indices(from_cc)[2:end]
    for i in 1:N
        # @info i
		mi = minimum(dismap[ind[i]])
		dis[i] = mi
		distmap[ind[i]] .= mi
		overlap[i] = sum(to_mask[ind[i]])
        # @inbounds dis[i] = minimum(dismap[ind[i]])
    end
    @assert(all(dis .>= 0))
	@debug "Distances $N"
    return dis, overlap, distmap
end

"""
	colocalize(xs, ys; metric="pearson", windowsize=3)

	Compute a colocalization metric `metric` on a window size of `windowsize^2`.

	Set `windowsize` to >= 3, odd. If set to -1, use the entire image, output will be single scalar.

	See `metrics` for metrics to use.
"""
function colocalize_2d(_xs, _ys; metric="pearson", windowsize=3)
	xs = Float64.(_xs)
	ys = Float64.(_ys)
	@info "Coloc with window $windowsize  and metric $metric for input $(size(xs))"
    X, Y = size(xs)
	if size(xs)!=size(ys)
		@error "Dimensions $(size(xs)) $(size(ys))"
		throw(ArgumentError("Dimensions $(size(xs)) $(size(ys))"))
	end
	if windowsize == -1
		return colocalize_nowindow(_xs, _ys; metric=metric)
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

function colocalize(_xs::Array{T, 3}, _ys::Array{T, 3}; metric="pearson", windowsize=3) where T <: Number
	@info "3D"
	xs = Float64.(_xs)
	ys = Float64.(_ys)
	@info "Coloc with window $windowsize  and metric $metric for input $(size(xs))"
    X, Y, Z = size(xs)
	if size(xs)!=size(ys)
		@error "Dimensions $(size(xs)) $(size(ys))"
		throw(ArgumentError("Dimensions $(size(xs)) $(size(ys))"))
	end
	if windowsize == -1
		return colocalize_nowindow(_xs, _ys; metric=metric)
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
    mf = metrics[metric]
    result = aszero(xs)
	@showprogress for zi in k+1:(Z-k)
		for yi in k+1:(Y-k)
			for xi in k+1:(X-k)
    			result[xi, yi, zi] = nanz(mf(xs[xi-k:xi+k, yi-k:yi+k, zi-k:zi+k], ys[xi-k:xi+k, yi-k:yi+k, zi-k:zi+k]))
			end
    	end
	end
    return result
end

function colocalize(_xs::Array{T, 2}, _ys::Array{T, 2}; metric="pearson", windowsize=3) where T <: Number
	@info "2D"
	return colocalize_2d(_xs, _ys; metric=metric, windowsize=windowsize)
end

function colocalize(_xs::Matrix, _ys::Matrix; metric="pearson", windowsize=3)
	@info "2D"
	return colocalize_2d(_xs, _ys; metric=metric, windowsize=windowsize)
end

function colocalize(_xs, _ys; metric="pearson", windowsize=3)
	@error "NOT IMPLEMENTED"
	throw(ArgumentError("Invalid arugment type $(type(_xs))"))
end

function colocalize_nowindow(_xs, _ys; metric="pearson")
	mf = metrics[metric]
	return nanz(mf(Float64.(_xs),Float64.(_ys)))
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

metrics = Dict([("pearson", coloc_pearson), ("spearman", coloc_spearman), ("jaccard", coloc_jaccard), ("manders", coloc_manders), ("m1", coloc_m1), ("sorensen", coloc_sorensen), ("m2", coloc_m2), ("haussdorff_max", haussdorff_max), ("haussdorff_mean", haussdorff_mean)])

"""
	list_metrics()

	Returns a list of the supported metrics
"""
function list_metrics()
	return keys(metrics) |> collect
end

"""
	metrics_iterator()

	Returns an iterator over the supported metrics, in String / Function pairs.

	Metric functions expect equal shaped array, and return a scalar value
"""
function metrics_iterator()
	return pairs(copy(metrics))
end

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
