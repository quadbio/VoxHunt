context('Tests for all basic VoxHunt functions')

library(testthat)
library(voxhunt)

#### ABA Voxel Maps ####
test_that('load_aba_data warns if files are not present.', {
    expect_warning(load_aba_data('/'))
    expect_warning(load_aba_data('../../../voxhunt_rds/'), NA)
})

test_that('data can be loaded.', {
    expect_error(data('example_seurat'), NA)
    expect_error(data('voxel_meta'), NA)
    expect_error(data('brainspan'), NA)
})

data('example_seurat')
load_aba_data('../../../voxhunt_rds/')

test_that('expression plots do not throw errors.', {
    expect_error(plot_expression('E13', gene = c('NEUROD6', 'DLX2')), NA)
    expect_error(plot_expression('E13', gene = c('NEUROD6', 'DLX2'), view = 'slice'), NA)
    expect_error(plot_expression_3d('E15', gene = c('NEUROD6', 'DLX2')), NA)
})

test_that('annotation plots do not throw errors.', {
    expect_error(plot_annotation('E13'), NA)
    expect_error(plot_annotation('E18'), NA)
})

test_that('voxel_map does not throw errors.', {
    expect_error(voxel_map(example_seurat, 'E13', genes_use = c('NEUROD6', 'DLX2')), NA)
})

vm <- voxel_map(example_seurat, 'E13', genes_use = c('NEUROD6', 'DLX2', 'LHX2'))

test_that('plot_map does not throw errors.', {
    expect_error(plot_map(vm), NA)
    expect_error(plot_map(vm, view='slice'), NA)
})

test_that('plot_map_3d does not throw errors.', {
    expect_error(plot_map_3d(vm), NA)
    expect_error(plot_map_3d(vm, annotation_level='custom_2'), NA)
})

test_that('plot_structure_similarity does not throw errors.', {
    expect_error(plot_structure_similarity(vm, groups=example_seurat$cluster), NA)
})


