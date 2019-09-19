#' Map correlation matrix to brain annotations and plot the result
#'
#' @import dplyr
#'
#' @export
#'
project_to_slices <- function(corr_mat,
                              cell_groups = NULL,
                              slices = c(6, 11, 23, 28, 36),
                              annotation_level = 'level_1',
                              annotation_colors = many
                              ){

    cluster_cor <- aggregate_matrix(corr_mat, groups=cell_groups, fun=colMeans)

    plot_df <- cluster_cor %>%
        as.matrix() %>%
        as_tibble(rownames='voxel') %>%
        tidyr::gather(cluster, corr, -voxel)

    if (is.null(cell_groups)){
        grouping_levels <- ' '
    } else{
        grouping_levels <- levels(factor(cell_groups))
    }

    utils::data(voxel_meta, envir = environment())
    meta <- filter(voxel_meta, voxel%in%colnames(corr_mat)) %>%
        tibble::column_to_rownames('voxel')

    plot_df$struct <- meta[plot_df$voxel, ][[annotation_level]]
    plot_df$x <- meta[plot_df$voxel, ]$x
    plot_df$y <- meta[plot_df$voxel, ]$y
    plot_df$z <- meta[plot_df$voxel, ]$z

    plot_df <- plot_df %>%
        filter(x%in%slices) %>%
        mutate(slice=as.numeric(factor(x))) %>%
        mutate(cluster=factor(cluster, levels=grouping_levels))

    slice_plot(plot_df, struct_colors = annotation_colors, newpage = T)

    return(plot_df)
}
