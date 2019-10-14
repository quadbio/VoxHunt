#' Plot brain structure annotation
#'
#' @rdname plot_annotation
#' @export plot_annotation
#'
plot_annotation <- function(object = 'E13', ...){
    UseMethod(generic = 'plot_annotation', object = object)
}


#' Plot brain structure annotation
#'
#' @import dplyr
#'
#' @rdname plot_annotation
#' @export
#' @method plot_annotation default
#'
plot_annotation.default <- function(
    object = 'E13',
    annotation_level = 'custom_2',
    annotation_colors = many,
    slices = NULL,
    show_coordinates = F,
    show_legend = F,
    alpha = 0.5
){
    utils::data(voxel_meta, envir = environment())
    age <- stage_name(object)
    meta <- filter(voxel_meta, stage == age)
    if (!is.null(slices)){
        meta <- filter(meta, z %in% slices)
    }
    p <- annotation_plot(
        annot_df = meta,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        show_coordinates = show_coordinates,
        show_legend = show_legend,
        alpha = alpha
    )
    return(p)
}


#' Plot brain structure annotation
#'
#' @import dplyr
#'
#' @rdname plot_annotation
#' @export
#' @method plot_annotation VoxelMap
#'
plot_annotation.VoxelMap <- function(
    object,
    annotation_level = 'custom_2',
    annotation_colors = many,
    slices = NULL,
    show_coordinates = F,
    show_legend = F,
    alpha = 0.5
){
    meta <- object$voxel_meta
    if (!is.null(slices)){
        meta <- filter(meta, z %in% slices)
    }
    p <- annotation_plot(
        annot_df = meta,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        show_coordinates = show_coordinates,
        show_legend = show_legend,
        alpha = alpha
    )
    return(p)
}

#' Plot brain structure annotation
#'
#' @import ggplot2
#'
annotation_plot <- function(
    annot_df,
    annotation_level = 'custom_2',
    annotation_colors = many,
    show_coordinates = F,
    show_legend = F,
    alpha = 0.5
){
    age <- unique(annot_df$stage)
    p <- ggplot(annot_df, aes_string('x', 'y', fill = annotation_level)) +
        geom_tile(alpha = alpha) +
        scale_x_continuous(breaks = seq(0, 100, 2)) +
        scale_y_continuous(breaks = seq(-100, 100, 2)) +
        scale_fill_manual(values = annotation_colors) +
        labs(title = age)
    if (!show_coordinates){
        p <- p + theme_void()
    }
    if (!show_legend){
        p <- p + theme(legend.position = 'none')
    }
    return(p)
}




#' Plot voxel map of single cells to voxels
#'
#' @rdname plot_map
#' @export plot_map
#'
plot_map <- function(object, ...){
    UseMethod(generic = 'plot_map', object = object)
}


#' Plot voxel map of single cells to voxels
#'
#' @import ggplot2
#' @import dplyr
#'
#' @rdname plot_map
#' @export
#' @method plot_map VoxelMap
#'
plot_map.VoxelMap <- function(
    object,
    view = 'saggital',
    slices = c(6, 11, 23, 28, 36),
    groups = NULL,
    annotation_level = 'custom_2',
    annotation_colors = many,
    map_colors = gyrdpu,
    show_coordinates = F,
    show_legend = F,
    newpage = T
){
    possible_views <- c('saggital', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice')
    if (!view %in% possible_views){
        stop(cat(
            paste0('"', view, '" is not a valid view argument. Try one of these:\n'),
            paste(possible_views, collapse = ', ')
            ))
    }

    if (is.null(groups) & !is.null(object$cell_meta$group)){
        groups <- object$cell_meta$group
    }

    if (is.null(groups)){
        grouping_levels <- ' '
    } else{
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

    if (view == 'slice'){
        plot_df <- plot_df %>%
            filter(x%in%slices) %>%
            mutate(slice=factor(x))

        p <- slice_plot(
            slice_df = plot_df,
            annotation_colors = annotation_colors,
            map_colors = map_colors,
            newpage = newpage
        )
    } else {
        p <- mapping_plot(
            map_df = plot_df,
            slices = slices,
            view = view,
            map_colors = map_colors,
            newpage = T
        )
    }
    return(p)
}


#' Plot correlation to brain slices and structure annotation
#'
#' @import ggplot2
#' @import dplyr
#'
slice_plot <- function(
    slice_df,
    annotation_colors = many,
    map_colors = gyrdpu,
    newpage = T
){
    annot <- ggplot(slice_df, aes(z, y, fill = struct)) +
        geom_tile() +
        theme_void() +
        scale_fill_manual(values = annotation_colors) +
        facet_grid(slice~., switch = 'both') +
        theme(
            legend.position = 'none',
            strip.text = element_text(angle = 180, size = 10)
        )

    plots <- purrr::map(levels(factor(slice_df$group)), function(g){
        plot_df <- slice_df %>%
            filter(group==g) %>%
            dplyr::group_by(z, y, slice) %>%
            dplyr::summarize(corr = mean(corr, na.rm = T))

        p <- ggplot(plot_df, aes(-z, y, fill = corr, alpha = corr)) +
            geom_tile() +
            theme_void() +
            scale_fill_gradientn(colors = map_colors) +
            scale_alpha_continuous(range = c(0.5, 1)) +
            facet_grid(slice ~ .) +
            theme(
                legend.position = 'none',
                strip.text = element_blank(),
                plot.title = element_text(angle = 60, hjust = 0.5, vjust = 0, size = 10)
                ) +
            labs(title = g)
        return(p)
    })

    plots <- c(list(annot), plots)
    p <- egg::ggarrange(plots = plots, nrow = 1, newpage = newpage)
    print(p)
}


#' Plot correlation to entire brain from saggital view
#'
#' @import ggplot2
#' @import dplyr
#'
mapping_plot <- function(
    map_df,
    view = 'saggital',
    slices = NULL,
    map_colors = rdpu,
    newpage = T
){
    plots <- map(levels(factor(map_df$group)), function(g){
        plot_df <- filter(map_df, group==g)
        if (view %in% c('coronal', 'x')){
            if (!is.null(slices)){
                plot_df <- filter(plot_df, x %in% slices)
            }
            p <- ggplot(plot_df, aes(z, y, fill = corr, alpha = corr))
        } else if (view %in% c('traverse', 'y')){
            if (!is.null(slices)){
                plot_df <- filter(plot_df, y %in% slices)
            }
            p <- ggplot(plot_df, aes(x, z, fill = corr, alpha = corr))
        } else {
            if (!is.null(slices)){
                plot_df <- filter(plot_df, z %in% slices)
            }
            p <- ggplot(plot_df, aes(x, y, fill = corr, alpha = corr))
        }
        p <- p +
            geom_tile() +
            theme_void() +
            theme(
                legend.position = 'none',
                plot.title = element_text(angle = 45, hjust = 0, vjust = 0, size = 8)
                ) +
            scale_fill_gradientn(colors = map_colors) +
            labs(title = g)
        return(p)
    })

    p <- egg::ggarrange(plots = plots, nrow = 1, newpage = newpage)
    print(p)
}


#' Plot gene expression across the mouse brain
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
plot_expression <- function(
    stage = 'E13',
    genes = NULL,
    view = 'saggital',
    slices = NULL
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
    }

    possible_views <- c('saggital', 'coronal', 'traverse', 'z' , 'x', 'y')
    if (!view %in% possible_views){
        stop(cat(
            paste0('"', view, '" is not a valid view argument. Try one of these:\n'),
            paste(possible_views, collapse = ', ')
        ))
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    voxel_meta <- DATA_LIST[[stage]]$row_meta

    if (view %in% c('coronal', 'x')){
        if (!is.null(slices)){
            voxel_meta <- filter(voxel_meta, x %in% slices)
        }
        f_plot <- function(x, g){
            ggplot(x, aes(z, y, fill=expr, alpha=expr)) +
                geom_tile(alpha=point_alpha) +
                dr_theme +
                labs(title=g) +
                feature_fill_scale +
                feature_theme
        }
    } else if (view %in% c('traverse', 'y')) {
        if (!is.null(slices)){
            voxel_meta <- filter(voxel_meta, y %in% slices)
        }
        f_plot <- function(x, g){
            ggplot(x, aes(x, z, fill=expr, alpha=expr)) +
                geom_tile(alpha=point_alpha) +
                dr_theme +
                labs(title=g) +
                feature_fill_scale +
                feature_theme
        }
    } else {
        if (!is.null(slices)){
            voxel_meta <- filter(voxel_meta, z %in% slices)
        }
        f_plot <- function(x, g){
            ggplot(x, aes(x, y, fill=expr, alpha=expr)) +
                geom_tile(alpha=point_alpha) +
                dr_theme +
                labs(title=g) +
                feature_fill_scale +
                feature_theme
        }
    }

    p <- feature_plot(
        expr_mat = voxel_mat,
        meta = voxel_meta,
        markers = genes,
        plot = f_plot
    )
    return(p)
}

#' Feature plot
#'
#' @import ggplot2
#' @import dplyr
#'
feature_plot <- function(
    expr_mat,
    meta,
    markers,
    plot = xy_plot,
    title = NULL,
    ncol = NULL,
    nrow = NULL,
    sort = T,
    scale = T
){
    meta$voxel <- as.character(meta$voxel)
    expr_markers <- get_markers(expr_mat, markers, scale=scale) %>%
        mutate(voxel=as.character(voxel)) %>%
        right_join(meta)
    plots <- map(unique(markers), function(g){
        x <- expr_markers %>%
            filter(gene==g)
        if (sort){
            x <- arrange(x, expr)
        }
        plot(x, g)
    })
    do.call(grid.arrange, args=c(plots, top=title, ncol=ncol, nrow=nrow))
}


