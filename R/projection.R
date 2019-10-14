#' Slice VoxelMap and plot result
#'
#' @rdname slice
#' @export slice
#'
slice <- function(object, ...){
    UseMethod(generic = 'slice', object = object)
}


#' Slice VoxelMap and plot result
#'
#' @import dplyr
#'
#' @rdname slice
#' @export
#' @method slice VoxelMap
#'
slice.VoxelMap <- function(
    object,
    slices = c(6, 11, 23, 28, 36),
    groups = NULL,
    annotation_level = 'level_1',
    annotation_colors = many
){
    if (is.null(groups) & !is.null(object$cell_meta$group)){
        groups <- object$cell_meta$group
    }

    if (is.null(groups)){
        grouping_levels <- ' '
    } else {
        grouping_levels <- levels(factor(groups))
    }

    cluster_cor <- aggregate_matrix(object$corr_mat, groups=groups, fun=colMeans)

    plot_df <- cluster_cor %>%
        as.matrix() %>%
        as_tibble(rownames='voxel') %>%
        tidyr::gather(group, corr, -voxel) %>%
        mutate(group=factor(group, levels=grouping_levels))

    meta <- column_to_rownames(object$voxel_meta, 'voxel')
    plot_df$struct <- meta[plot_df$voxel, ][[annotation_level]]
    plot_df$x <- meta[plot_df$voxel, ]$x
    plot_df$y <- meta[plot_df$voxel, ]$y
    plot_df$z <- meta[plot_df$voxel, ]$z

    plot_df <- plot_df %>%
        filter(x%in%slices) %>%
        mutate(slice=as.numeric(factor(x)))

    slice_plot(
        slice_df = plot_df,
        annotation_colors = annotation_colors,
        newpage = T
    )

    object$slices <- plot_df

    return(object)
}

