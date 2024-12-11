# Colocalization

A Julia package providing colocalization metrics for images and their sparse representations.

This package allows you to quickly run all metrics, and report the results both in image and CSV format.

Colocalization is used often in multichannel microscopy to quantify functional interaction between fluorescently marked proteins or subcellular organelles.
Note that colocalization in superresolution microscopy has to be very carefully applied, as with increasing precision no two objects can share the same location at the same time.

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/bencardoen/Colocalization.jl/tree/main.svg?style=svg&circle-token=50ed75938474a05f8c9ed7343d9d6134131f5519)](https://dl.circleci.com/status-badge/redirect/gh/bencardoen/Colocalization.jl/tree/main) [![codecov](https://codecov.io/gh/bencardoen/Colocalization.jl/branch/main/graph/badge.svg?token=50R4ZYYY1V)](https://codecov.io/gh/bencardoen/Colocalization.jl) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7552357.svg)](https://doi.org/10.5281/zenodo.7552357)

## Table of contents
1. [Installation](#installation)
2. [Usage](#usage)

    2.0 [Point Clouds](#pc)

    2.0 [Voxel based data](#vx)

        2.1 [Supported Metrics](#metrics)
        
        2.2 [Demo](#demo)

        2.3 [Documentation](#docs)

3. [Cite](#cite)
4. [FAQ](#faq)
5. [Related projects and tools](#related)
6. [Can you support metric X?](#support)


<a name="installation"></a>
## Installation
- Download [Julia](https://julialang.org/learning/getting-started/)
- [Julia + VSCode](https://blog.glcs.io/install-julia-and-vscode#heading-installing-a-specific-julia-version)
- Open a new VSCode window
- In the terminal, type `git clone https://github.com/bencardoen/Colocalization.jl.git` 
- change directory to `Colocalization.jl` which will now be a subdirectory
### Using as a package
Start Julia (in VSCode or Command line)
```bash
julia
```
In Julia
```julia
using Pkg
# Optionally, activate your environment
# Pkg.activate("path/to/your/environment")
Pkg.add(url="https://github.com/bencardoen/Colocalization.jl")
using Colocalization
```

### Cloning the repository
This assumes you have [Git](https://git-scm.com/downloads) installed and configured.
```bash
git clone https://github.com/bencardoen/Colocalization.jl.git
cd Colocalization.jl
julia --project=.
```
Then, in Julia:
```julia
using Pkg
Pkg.instantiate()
using Colocalization
```
That's it.

#### On Command line
Let's say you have 2 image files `a.tif` and `b.tif`.
```
julia --project=. scripts/colocalize.jl -f a.tif -s b.tif --outdir . --segment -w 3
```

<a name="usage"></a>
## Usage

<a name="pc"></a>
### Point cloud
```julia
julia scripts/colocalize_pointcloud.jl --first 1st.mat --second 2nd.mat --outdir "X"
```
This reads in SuperResNet files (MAT format) of two channels, 3D localizations.
The output will be a CSV file where each row describes the colocalization of **1** cluster in channel x to **5** clusters to channel y.
The columns are:
- channel: e.g. 1, 2
- channel name: the corresponding filename
- clusterid : this is the integer identifier used for this cluster in SRN
- centroid_{x,y,z} : the centroid location of this cluster
- distance_{1-5} : The distances to the nearest 5 objects in the other channel
- nearest_{1-5} : The cluster ids to the nearest 5 objects in the other channel
- channel_centroid_{x,y,z} : The centroid of **this** channel
- distance_to_centroid: the distance of this object's centroid to the channel's centroid (~ density/topology)
- radius: This is the radius of the circumscribed sphere of this cluster. If, for two clusters, you have R1 and R2 as radii, and their centroid to centroid distance is D12, then you can detect `overlap` = D12 < R1 + R2. This is an approximate measure, as the cluster can have weird shapes that skew the size of the circumscribed circle. 

#### SRN Specific data
You can extract specific fields from the SRN MAT file.
Start julia
```bash
julia --project=.
```
Then, for example to read the 3D points of the clusters, you can do

```julia
using Colocalization
data = load_SRN("file.mat")
points = data[4]
CSV.write("3D-points.csv", DataFrame(points, ["X", "Y", "Z"]))
```

<a name="vx"></a>
### Voxel 

<a name="metrics"></a>
### Supported Metrics 

You can get an up to date listing of the supported metrics by running the following code:
```julia
using Colocalization, Logging
@info list_metrics()
```
or access the actual functions:
```julia
for (name, metric) in metrics_iterator()
  @info name, metric
end
```

<a name="demo"></a>
### In silico example

Let's create 2 objects with variable levels of fluorescence labelling, that overlap by 50%.
```julia
using Images, Statistics, Distributions, Colocalization, ImageFiltering, Random
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
using SPECHT, ImageView
imshow(mosaicview( [SPECHT.tcolors([xs, ys]), SPECHT.tcolors([sx, sy]), SPECHT.tcolors([s2x, s2y])], nrow=1))
```

The visualzation snippet uses SPECHT and Imageview, if you don't have them:
```julia
using Pkg
Pkg.add("ImageView")
Pkg.add(url="https//github.com/bencardoen/SPECHT.jl")
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

<a name="docs"></a>
### Documentation
The documentation of the functions describes proper usage and meaning of parameters, to access it:
```julia
using Colocalization
?colocalize_all
```
The `?` key invokes Julia documentation, tools/IDES such as VSCode/Atom would have built in documentation panes.

 <a name="cite"></a>
## Cite
If you find this useful, consider citing:
```bibtext
@software{ben_cardoen_2023_7552357,
  author       = {Ben Cardoen},
  title        = {Colocalization.jl},
  month        = jan,
  year         = 2023,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.7552357},
  url          = {https://doi.org/10.5281/zenodo.7552357}
}
```

**Note** For the individual metrics, please cite the introducing author!!!.

<a name="faq"></a>
## FAQ
- To display the images, you need to install ImageView
```julia
using Pkg
Pkg.add("ImageView")
```
If you have any problems or suggestions, please create an [issue](https://github.com/bencardoen/Colocalization.jl/issues/new/choose)

 <a name="related"></a>
## Related software
FiJi:
- [https://imagej.net/plugins/coloc-2](https://imagej.net/plugins/coloc-2)
- [https://github.com/fiji/Colocalisation_Analysis](https://github.com/fiji/Colocalisation_Analysis)

This package would not be possible without the [Julia Images ecosystem](https://juliaimages.org/latest/)

<a name="support"></a>
## Can you support Metric X?
Sure, please create an [issue](https://github.com/bencardoen/Colocalization.jl/issues/new/choose) describing the metric mathematically, ideally accompanied by the introducing paper.
