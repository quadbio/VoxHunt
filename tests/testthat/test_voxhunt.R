library(testthat)
library(voxhunt)

data('example_seurat')

test_that('load_aba_data warns if files are not present.', {
    expect_warning(load_aba_data('/'))
    expect_error(load_aba_data('../voxhunt_rds/'), NA)
})





