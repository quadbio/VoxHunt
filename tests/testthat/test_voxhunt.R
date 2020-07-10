context('Tests for all basic VoxHunt functions')

library(testthat)
library(voxhunt)

data('example_seurat')

test_that('load_aba_data warns if files are not present.', {
    expect_warning(load_aba_data('/'))
    expect_warning(load_aba_data('../../../voxhunt_rds/'), NA)
})

test_that('data can be loaded.', {
    expect_error(data('example_seurat'), NA)
    expect_error(data('voxel_meta'), NA)
    expect_error(data('brainspan'), NA)
})

load_aba_data('../../../voxhunt_rds/')
data('example_seurat')

test_that('expression plots do not throw errors.', {
    expect_error(plot_expression('E13', gene = c('NEUROD6', 'DLX2')), NA)
    expect_error(plot_expression('E13', gene = c('NEUROD6', 'DLX2'), view = 'slice'), NA)
    expect_error(plot_expression_3d('E13', gene = c('NEUROD6', 'DLX2')), NA)
})

test_that('voxel_map does not throw errors.', {
    expect_error(voxel_map(example_seurat, 'E13', genes_use = c('NEUROD6', 'DLX2')), NA)
})

vm <- voxel_map(example_seurat, 'E13', genes_use = c('NEUROD6', 'DLX2'))

test_that('plot_map does not throw errors.', {
    expect_error(plot_map(vm), NA)
    expect_error(plot_map(vm, view='slice'), NA)
})


