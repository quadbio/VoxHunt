# VoxHunt <img src="figures/logo.png" align="right" width="120" />


You want to validate your new brain organoid protocol, find out which cell types emerge in your current batch of cerebral organoids or just find a very specific markers for a tiny brain structure? 

You have come to the right place. 


## Introduction 

Brain organoids are complex and can contain cells at various stages of differentiation from different brain structures. Single cell genomic methods provide powerful approaches to explore cell composition, differentiation trajectories, gene regulation, and genetic perturbations in brain organoid systems. VoxHunt is a handy little tool to assess brain organoid patterning, developmental state, and cell composition through systematic comparisons to three-dimensional in situ hybridization data from the Allen Brain Atlas.

<img src="figures/abstract.png" align="center" />


## Installation

Presto, one of VoxHunt's dependencies is not on CRAN and has to be installed from GitHub:
```{r}
# install.packages('devtools')
devtools::install_github('immunogenomics/presto')
```
Once Presto is installed, you can install VoxHunt with
```{r}
devtools::install_github('quadbiolab/voxhunt')
```

## Quick start

If you have a `seurat_object` with single cell transcriptomic data of your organoid ready, you can start right away with projecting them to the brain:
```{r}
library(voxhunt)
genes_use <- voxhunt::variable_genes('E13', 300)$gene
vox_map <- voxhunt::voxel_map(seurat_object, genes_use=genes_use)
voxhunt::plot_map(vox_map)
voxhunt::plot_map(vox_map, view='slice')
```
Here, we select the 300 most variable features from the E13.5 mouse brain ISH data and use them to calculate similarity maps for your organoid cells. We then plot these maps in the saggital view and as coronal sections. 

If you want to find out about VoxHunt's functionality, have a look into our [vignette](http://htmlpreview.github.io/?https://github.com/quadbiolab/VoxHunt/blob/master/vignettes/getting_started.html).





