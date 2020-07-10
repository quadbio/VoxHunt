library(testthat)
library(voxhunt)

data('example_seurat')
load_aba_data('../voxhunt_rds/')
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


