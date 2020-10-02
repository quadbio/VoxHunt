library(testthat)
library(tidyverse)
library(patchwork)
library(voxhunt)
library(Seurat)

data('example_seurat')
load_aba_data('../voxhunt_data/')
plot_expression('E13', 'NEUROD6')
plot_expression('E15', 'NEUROD6', view='sagittal')


markers <- structure_markers('E13')
top10 <- markers %>%
    group_by(group) %>%
    top_n(10, auc) %>%
    {unique(.$gene)}


pb_voxmap <- voxel_map(example_seurat, 'E13', genes_use = top10, pseudobulk_groups = T, group_name = 'cluster')
plot_map(pb_voxmap)

voxmap <- voxel_map(example_seurat, 'E13', genes_use = top10)
plot_map(voxmap, groups = example_seurat$cluster)
plot_map(voxmap, groups = example_seurat$cluster, view='slice')
plot_structure_similarity(voxmap, groups = example_seurat$cluster)

### BrainSpan
pb_bsmap <- brainspan_map(example_seurat, genes_use = top10, pseudobulk_groups = T, group_name = 'cluster')
plot_map(pb_bsmap)

bsmap <- brainspan_map(
    example_seurat,
    genes_use = top10,
    group_name = 'cluster',
    pseudobulk_groups = F
)
plot_map(bsmap)
plot_structure_similarity(bsmap)

### La Manno Mousebrain
load_mousebrain_data('/Volumes/treutlein/PUBLIC_DATA/published/single_cell/2020_linnarson_dev_mouse_brain/')
pb_mbmap <- mousebrain_map(example_seurat, genes_use = top10, pseudobulk_groups = T, group_name = 'cluster')
plot_structure_similarity(pb_mbmap, cluster=T)
plot_structure_similarity(pb_mbmap, cluster=F)
plot_map(pb_mbmap)

mbmap <- mousebrain_map(example_seurat, genes_use = top10)
plot_structure_similarity(mbmap, groups = example_seurat$cluster, cluster=F)
plot_map(mbmap, groups = example_seurat$cluster)

### Deconvolution
data('example_pseudobulk')

markers <- structure_markers('E13', annotation_level = 'custom_2')

involve_regions <- c(
    'pallium', 'subpallium',
    'diencephalon', 'midbrain', 'hypothalamus',
    'prepontine hindbrain', 'pontine hindbrain', 'medullary hindbrain'
)

top15 <- markers %>%
    filter(group%in%involve_regions) %>%
    group_by(group) %>%
    top_n(15, auc) %>%
    {unique(.$gene)}

top50 <- markers %>%
    filter(group%in%involve_regions) %>%
    group_by(group) %>%
    top_n(50, auc) %>%
    {unique(.$gene)}

prop_df <- deconvolute(
    example_pseudobulk[top50, ], top15,
    involve_regions = involve_regions,
    pseudo_tpm = T
)

ggplot(prop_df, aes(sample, prop, fill=factor(struct, levels=struct_names_custom2))) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=struct_colors_custom2) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0))






