library(testthat)
library(tidyverse)
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

bsmap <- brainspan_map(example_seurat, genes_use = top10)
plot_map(bsmap, groups = example_seurat$cluster)

load_mousebrain_data('/Volumes/treutlein/PUBLIC_DATA/published/single_cell/2020_linnarson_dev_mouse_brain/')
mbmap <- mousebrain_map(example_seurat, genes_use = top10)


#### Reformat mouse conversion ###
m2h <- read_tsv('/Volumes/jfleck/resources/mappings/human2mouse.tsv')
conversion <- m2h
colnames(conversion)[1:5] <-c('homology_id', 'organism', 'ncbi_taxom_id', 'symbol', 'entrez_id')
conversion <- conversion %>%
    group_by(homology_id) %>%
    filter(length(unique(organism))==2)
conversion_h <- conversion %>%
    select(organism, symbol) %>%
    filter(organism=='human') %>%
    select('human_symbol'=symbol)
conversion_m <- conversion %>%
    select(organism, symbol) %>%
    filter(organism!='human') %>%
    select('mouse_symbol'=symbol)
conversion_both <- inner_join(conversion_h, conversion_m)
