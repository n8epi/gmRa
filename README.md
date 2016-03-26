# gmRa
Geometric Multiresolution Analysis in R

## Synopsis
Geometric Multiresolution Analysis [(Allard, Chen, Maggioni, 2012)](http://arxiv.org/pdf/1105.4924.pdf) is a data science technique for constructing efficient representations of high-dimensional datasets. Essentially, the construction has two stages:

1. Cluster, or partition the data into similar groups
2. Fit low-dimensional hyperplanes to each data cluster/partition

In this implementation, [cover trees](http://hunch.net/~jl/projects/cover_tree/cover_tree.html) are used to produce a subset of the data that "covers" the dataset at a particular scale and the data is partitioned using the [Voronoi regions](https://en.wikipedia.org/wiki/Voronoi_diagram) produced by this subset. Fitting of hyperplanes in each Voronoi region is carried out by first computing the local mean and covariance matrix, and identifying the unique hyperplane of dimension d containing the mean and parallel to the top d eigenvectors of the covariance matrix. Strong theoretical justifications for this procedure can be found in [Maggioni, Minsker, and Strawn, 2016](http://www.jmlr.org/papers/v17/maggioni16a.html).

This code is not fully optimized, but should be useful for educational and experimental purposes.

## License
gmRa is released under the MIT license

