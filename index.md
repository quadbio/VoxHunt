# VoxHunt <img src="man/figures/logo.png" align="right" width="150" />

VoxHunt is a package for assessing brain organoid patterning, developmental state, and cell composition through systematic comparisons of single cell transcriptomes to three-dimensional in situ hybridization data from the Allen Brain Atlas and a number of other useful reference datasets.


## Features

VoxHunt currently implements the following features:
* Exploration of gene expression in the developing mouse brain
* Discovery of brain structure markers
* Comparison of scRNA-seq data to
    * ISH data from the [Allen Developing Mouse Brain Atlas](https://developingmouse.brain-map.org/) 
    * bulk RNA-seq data of brain regions in the developing human brain ([BrainSpan](https://www.brainspan.org/))
    * scRNA-seq data from the developing mouse brain ([La Manno & Siletti et al.](https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1))
* Visualization of spatial similarity patterns
* Annotation of cells with brain regional identities
* Deconvolution of bulk RNA-seq data based on brain structure expression patterns

<a/>
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

## Citation

If you find VoxHunt useful for your research, please consider citing our paper:

[Fleck, JS, Sanchís-Calleja F et al., Resolving organoid brain region identities by mapping single-cell genomic data to reference atlases. Cell Stem Cell, 2020, https://doi.org/10.1016/j.stem.2021.02.015](https://doi.org/10.1016/j.stem.2021.02.015)

```
@Article{,
    title = {Resolving organoid brain region identities by mapping single-cell genomic data to reference atlases},
    author = {Jonas Simon Fleck and Fátima Sanchís-Calleja and Zhisong He and Margozata Santel and Michael James Boyle and J. Gray Camp and Barbara Treutlein},
    journal = {Cell Stem Cell},
    year = {2021},
    url = {https://doi.org/10.1016/j.stem.2021.02.015},
    doi = {10.1016/j.stem.2021.02.015}
}
```

