library(testthat)
library(voxhunt)

data('example_seurat')
load_aba_data('../voxhunt_rds/')
plot_expression('E13', 'NEUROD6')
plot_expression('E15', 'NEUROD6', view='sagittal')



