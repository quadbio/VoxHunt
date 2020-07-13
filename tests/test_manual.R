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

data('human2mouse')


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
human2mouse <- inner_join(conversion_h, conversion_m)
save(human2mouse, file='data/human2mouse.RData')


dir <- '/Volumes/treutlein/PUBLIC_DATA/published/single_cell/2020_linnarson_dev_mouse_brain/'

all_path <- file.path(dir, 'dev_all.loom')
all_loom <- hdf5r::H5File$new(all_path)
agg_path <- file.path(dir, 'dev_all.agg.loom')
agg_loom <- hdf5r::H5File$new(agg_path)

all_tsne <- all_loom[['col_attrs/TSNE']][,]
all_region <- all_loom[['col_attrs/Region']][]
all_class <- all_loom[['col_attrs/Class']][]
all_cluster <- all_loom[['col_attrs/Clusters']][]

all_meta <- tibble::tibble(
    tSNE1 = all_tsne[1, ],
    tSNE2 = all_tsne[2, ],
    region = all_region,
    class = all_class,
    cluster = as.character(all_cluster)
)

agg_expression <- agg_loom[['matrix']][,]
colnames(agg_expression) <- make.names(stringr::str_to_upper(agg_loom[['row_attrs/Gene']][]))
rownames(agg_expression) <- agg_loom[['col_attrs/Clusters']][]

MOUSEBRAIN_DATA <<- list(
    matrix = Matrix(agg_expression, sparse=T),
    meta = all_meta
)

all_loom$close_all()
agg_loom$close_all()


