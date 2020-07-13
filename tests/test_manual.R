library(testthat)
library(tidyverse)
library(voxhunt)

data('example_seurat')
load_aba_data('../voxhunt_data/')
plot_expression('E13', 'NEUROD6')
plot_expression('E15', 'NEUROD6', view='sagittal')


markers <- structure_markers('E13')
top10 <- markers %>%
    group_by(group) %>%
    top_n(10, auc) %>%
    {unique(.$gene)}

voxmap <- voxel_map(example_seurat, 'E13', genes_use = top10)
plot_map(voxmap, groups = example_seurat$cluster)
plot_map(voxmap, groups = example_seurat$cluster, view='slice')
plot_structure_similarity(voxmap, groups = example_seurat$cluster)

bsmap <- brainspan_map(example_seurat, genes_use = top10)
plot_map(bsmap, groups = example_seurat$cluster)

load_mousebrain_data('/Volumes/treutlein/PUBLIC_DATA/published/single_cell/2020_linnarson_dev_mouse_brain/')
mbmap <- mousebrain_map(example_seurat, genes_use = top10)

plot_structure_similarity(mbmap, groups = example_seurat$cluster, cluster=F)

plot_map(mbmap, groups = example_seurat$cluster)


