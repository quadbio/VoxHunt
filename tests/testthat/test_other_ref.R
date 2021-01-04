context('Tests for other references.')

library(testthat)
library(voxhunt)

data('example_seurat')

#### BrainSpan ####
test_that('brainspan map does not throw errors.', {
    expect_error(brainspan_map(example_seurat, 'E13', genes_use = c('NEUROD6', 'DLX2')), NA)
})

bsmap <- brainspan_map(example_seurat, genes_use = c('NEUROD6', 'DLX2'))

test_that('brainspan maps can be plotted.', {
    expect_error(plot_map(bsmap), NA)
    expect_error(plot_structure_similarity(bsmap), NA)
})


#### Mousebrain ####
test_that('brainspan map does not throw errors.', {
    expect_error(mousebrain_map(example_seurat, genes_use = c('NEUROD6', 'DLX2'), pseudobulk_groups = T, group_name = 'cluster'), NA)
    expect_error(mousebrain_map(example_seurat, genes_use = c('NEUROD6', 'DLX2')), NA)
})

mbmap <- mousebrain_map(example_seurat, genes_use = c('NEUROD6', 'DLX2'), pseudobulk_groups = T, group_name = 'cluster')

test_that('brainspan maps can be plotted.', {
    expect_error(plot_map(mbmap), NA)
    expect_error(plot_structure_similarity(mbmap), NA)
})

