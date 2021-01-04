context('Tests for deconvolution')

library(testthat)
library(voxhunt)

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
    filter(gene%in%rownames(example_pseudobulk)) %>%
    group_by(group) %>%
    top_n(15, auc) %>%
    {unique(.$gene)}

top50 <- markers %>%
    filter(group%in%involve_regions) %>%
    filter(gene%in%rownames(example_pseudobulk)) %>%
    group_by(group) %>%
    top_n(50, auc) %>%
    {unique(.$gene)}

test_that('deconvolute does not throw errors.', {
    expect_error(deconvolute(
        example_pseudobulk[top50, ], top15,
        involve_regions = involve_regions,
        pseudo_tpm = T
    ), NA)
})

test_that('deconvolute returns expected values.', {
    prop_df <- deconvolute(
        example_pseudobulk[top50, ], top15,
        involve_regions = involve_regions,
        pseudo_tpm = T
    )
    expect_gte(prop_df[prop_df$sample=='kucg_HipSci_1' & prop_df$struct=='pallium', 'prop'], 0.8)
    expect_gte(prop_df[prop_df$sample=='kucg_HipSci_2' & prop_df$struct=='pallium', 'prop'], 0.9)
    expect_gte(prop_df[prop_df$sample=='h409B2_67d_org1' & prop_df$struct=='diencephalon', 'prop'], 0.2)
    expect_gte(prop_df[prop_df$sample=='wibj_64d_org1' & prop_df$struct=='subpallium', 'prop'], 0.4)
})

