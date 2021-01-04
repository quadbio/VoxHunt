![Build](https://github.com/quadbiolab/voxhunt/workflows/build/badge.svg)

# VoxHunt <img src="man/figures/logo.png" align="right" width="120" />

You want to validate your new brain organoid protocol, find out which cell types emerge in your current batch of cerebral organoids or just find a very specific marker for a tiny brain structure? VoxHunt might be what you are looking for.


## Introduction

Brain organoids are complex and can contain cells at various stages of differentiation from different brain structures. Single cell genomic methods provide powerful approaches to explore cell composition, differentiation trajectories, gene regulation, and genetic perturbations in brain organoid systems. VoxHunt is a handy tool to assess brain organoid patterning, developmental state, and cell composition through systematic comparisons of single cell transcriptomes to three-dimensional in situ hybridization data from the Allen Brain Atlas.

<img src="man/figures/abstract.png" align="center" />


## Installation

[Presto](https://github.com/immunogenomics/presto), one of VoxHunt's dependencies is not on CRAN and has to be installed from GitHub:

```r
# install.packages('devtools')
devtools::install_github('immunogenomics/presto')
```

If you want to perform deconvolution, you also have to install EPIC manually:

```r
devtools::install_github('GfellerLab/EPIC')
```

Once all external dependencies are installed, you can install VoxHunt with
```r
devtools::install_github('quadbiolab/voxhunt')
```

In addition to the R package itself, you'll also need to download the the ABA expression data from [here](http://dx.doi.org/10.17632/g4xg38mwcn.2).

## Quick start

If you have a `seurat_object` with single cell transcriptomic data of your organoid ready, you can start right away with projecting them to the brain:

```r
# Load VoxHunt
library(voxhunt)

# Point VoxHunt to ABA expression data
load_aba_data('~/path/to/data')

# Find 300 most variable genes from the E13.5 mouse brain
genes_use <- variable_genes('E13', 300)$gene

# Calculate the similarity map of a seurat object to the E13.5 mouse brain
vox_map <- voxel_map(seurat_object, genes_use=genes_use)

# Plot the result
plot_map(vox_map)
```
After loading VoxHunt, we first point it to the directory with the unpacked [ABA expression data](http://dx.doi.org/10.17632/g4xg38mwcn.1). Then, we select the 300 most variable features from the E13.5 mouse brain ISH data and use them to calculate similarity maps for your organoid cells. Finally, we plot these maps in the sagittal view.

## Mapping to other references

In addition to the ABA, we have also included functions that facilitate mapping to the [BrainSpan data](https://www.brainspan.org/) (human microdissected brain) as well as the recently released developing mouse brain scRNA-seq atlas ([La Manno & Siletti et al.](https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1)).

## More

If you want to find out more about VoxHunt's functionality, have a look at vignettes on our [website](https://quadbiolab.github.io/VoxHunt/index.html).

## Citation

If you find VoxHunt useful for your research, please consider citing our preprint:

[Fleck, J.S., et al., Resolving brain organoid heterogeneity by mapping single cell genomic data to a spatial reference. bioRxiv, 2020, https://doi.org/10.1101/2020.01.06.896282 ](https://www.biorxiv.org/content/10.1101/2020.01.06.896282v1)
