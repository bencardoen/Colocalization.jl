# Colocalization

Provides colocalization metrics for images and their sparse representations.

The metrics used are (for now) standard metrics, which I needed for a different project as comparison/baseline. 
I'll update as I come across more metrics.

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/bencardoen/Colocalization.jl/tree/main.svg?style=svg&circle-token=50ed75938474a05f8c9ed7343d9d6134131f5519)](https://dl.circleci.com/status-badge/redirect/gh/bencardoen/Colocalization.jl/tree/main) [![codecov](https://codecov.io/gh/bencardoen/Colocalization.jl/branch/main/graph/badge.svg?token=50R4ZYYY1V)](https://codecov.io/gh/bencardoen/Colocalization.jl)


## In silico example
spearman", "m2", "m1", "jaccard", "manders", "sorensen", "pearson"]
Let's create 2 objects with variable levels of fluorescence labelling, that overlap by 50%.
```julia
using ImageView, Images, Statistics, Distributions, SPECHT, Colocalization, ImageFiltering, Random
X, Y = 100, 100
xs = zeros(X, Y)
ys = zeros(X, Y)
xs[40:50, 40:50] .= rand(11, 11)
ys[45:55, 45:55] .= rand(11, 11)
sx = ImageFiltering.imfilter(xs, ImageFiltering.Kernel.gaussian((3, 3)))
sy = ImageFiltering.imfilter(ys, ImageFiltering.Kernel.gaussian((3, 3)))
```
We'll add some noise to make things realistic
```julia
s2x = copy(sx)
s2y = copy(sy)
s2x .+= rand(100, 100) ./ 10
s2y .+= rand(100, 100) ./ 10
```
View the results
```julia
imshow(mosaicview( [SPECHT.tcolors([xs, ys]), SPECHT.tcolors([sx, sy]), SPECHT.tcolors([s2x, s2y])], nrow=1))
```
This should produce something like this image

![demo.png](demo.png)

Now, we compute all coloc metrics
```julia
results = colocalize_all(s2x, s2y)
```
Let's view the results, the metrics from left to right are: `spearman, m2, m1, jaccard, manders, sorensen, pearson`
```
mv = mosaicview([abs.(results[k]) for k in keys(results)], nrow=1)
imshow(mv)
```

![demo.png](resultsnseg.png)

Clearly, the noise is throwing a wrench in things. Metrics like Jacard, M1 and so forth expect segmented images to work on.
Let's do a quick segmentation.
```julia
xt = otsu_threshold(s2x)
yt = otsu_threshold(s2y)
s2x[s2x.<xt] .= 0
s2y[s2y.<yt] .= 1
results = colocalize_all(s2x, s2y)
mv = mosaicview([abs.(results[k]) for k in keys(results)], nrow=1)
imshow(mv)
```
Which should produce something like the below image.
![demo.png](resultseg.png)

You can pass parameters, such as windowsize (default 3x3) to the function, please see their docstrings:
```julia
using Colocalization
?colocalize_all
```

## Installation
- Download [Julia](https://julialang.org/learning/getting-started/)

### Cloning the repository
```bash
git clone https://github.com/bencardoen/Colocalization.jl.git
cd Colocalization.jl
julia --project=.
```
in Julia
```julia
using Pkg
Pkg.instantiate()
```
That's it.

### Using as a package
```bash
julia
```
In Julia
```julia
using Pkg
Pkg.add(url="https://github.com/bencardoen/ERGO.jl")
Pkg.add(url="https://github.com/bencardoen/SPECHT.jl")
Pkg.add(url="https://github.com/bencardoen/Colocalization.jl")
using Colocalization
```
