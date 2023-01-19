# Colocalization

Provides colocalization metrics for images and their sparse representations.

## Installation

## Usage
```julia
using Colocalization
using Random
using Images
using ImageView
X = 100
Y = 200
# Create 2 random images
xs = Random.rand(X, Y)
ys = Random.rand(X, Y) .+ .2
res = colocalize_all(xs, ys)
# res is a dictionary with each metric being the key to the results
for metric in keys(Colocalization.metrics)
        coloc_m = res[metric]
end
# Let's summarize the results in a dataframe
df = summarize_colocalization(res, "f1", "f2")
```
Colocalization is defined on a window of 3x3 by default, you can change this:
```julia
res = colocalize_all(xs, ys; windowsize=11)
```

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
