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
        scale_fill_manual(values = annotation_colors, na.value='gray') +
        labs(title = age, fill= 'Structure') +
        theme_bw() +
        dr_guides
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
    view = 'sagittal',
    slices = NULL,
    groups = NULL,
    annotation_level = 'custom_2',
    annotation_colors = many,
    show_group = NULL,
    map_colors = gyrdpu_flat,
    show_coordinates = F,
    show_legend = F
){
    possible_views <- c('sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D')
    if (!view %in% possible_views){
        stop(cat(
            paste0('"', view, '" is not a valid view argument. Try one of these:\n'),
            paste(possible_views, collapse = ', '), '\n'
            ))
    }

    plot_df <- summarise_groups(object, groups)

    if (view == 'slice'){
        if (is.null(slices)){
            stop('To view slices you have to define a set of sections with the "slices" argument.')
        }
        plot_df <- plot_df %>%
            filter(x%in%slices) %>%
            mutate(slice=factor(x))

        p <- slice_plot(
            slice_df = plot_df,
            annotation_colors = annotation_colors,
            annotation_level = annotation_level,
            map_colors = map_colors
        )
    } else {
        p <- mapping_plot(
            map_df = plot_df,
            slices = slices,
            view = view,
            map_colors = map_colors
        )
    }
    return(p)
}

#' Plot voxel map of single cells to voxels in 3D
#'
#' @rdname plot_map_3d
#' @export plot_map_3d
#'
plot_map_3d <- function(object, ...){
    UseMethod(generic = 'plot_map_3d', object = object)
}


#' Plot voxel map of single cells to voxels in 3D
#'
#' @import ggplot2
#' @import dplyr
#'
#' @rdname plot_map_3d
#' @export
#' @method plot_map_3d VoxelMap
#'
plot_map_3d.VoxelMap <- function(
    object,
    groups = NULL,
    show_group = NULL,
    annotation_level = NULL,
    annotation_colors = many,
    map_colors = inferno,
    sizes = c(10, 1000)
){
    possible_views <- c('sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice')
    if (!view %in% possible_views){
        stop(cat(
            paste0('"', view, '" is not a valid view argument. Try one of these:\n'),
            paste(possible_views, collapse = ', '), '\n'
        ))
    }

    plot_df <- summarise_groups(object, groups)

    if (is.null(show_group) & is.null(groups)){
        show_group <- levels(factor(groups))[1]
    } else if (is.null(show_group) & !is.null(groups)){
        show_group <- levels(factor(object$cell_meta$group))[1]
    }

    plot_df <- filter(plot_df, group==show_group) %>%
        mutate(intensity=corr)
    plot_df_3d <- plot_df %>%
        mutate(z=-z+max(z)*2) %>%
        bind_rows(plot_df, .id = 'hemisphere')
    p <- three_dim_plot(
        int_df = plot_df_3d,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        intensity_colors = map_colors,
        sizes = sizes
    )
    return(p)
}


#' Plot correlation to brain slices and structure annotation
#'
#' @import ggplot2
#' @import dplyr
#'
slice_plot <- function(
    slice_df,
    annotation_level = 'custom_2',
    annotation_colors = many,
    map_colors = gyrdpu
){
    annot <- ggplot(slice_df, aes_string('z', 'y', fill = annotation_level)) +
        geom_tile() +
        theme_void() +
        scale_fill_manual(values = annotation_colors) +
        facet_grid(slice~., switch = 'both') +
        theme_bw() +
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
            theme_bw() +
            theme(
                legend.position = 'none',
                strip.text = element_blank(),
                plot.title = element_text(angle = 60, hjust = 0, vjust = 0, size = 10)
                ) +
            labs(title = g)
        return(p)
    })
    plots <- c(list(annot), plots)
    p <- egg::ggarrange(plots = plots, nrow = 1, draw = F)
    return(p)
}


#' Plot correlation to entire brain from sagital view
#'
#' @import ggplot2
#' @import dplyr
#'
mapping_plot <- function(
    map_df,
    view = 'sagittal',
    slices = NULL,
    map_colors = rdpu
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
            theme_bw() +
            theme_void() +
            theme(
                legend.position = 'none',
                plot.title = element_text(size = 8)
                ) +
            scale_fill_gradientn(colors = map_colors) +
            labs(title = g)
        return(p)
    })
    p <- egg::ggarrange(plots = plots, nrow = 1, draw = F)
    return(p)
}


#' Plot gene expression across the mouse brain
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
plot_expression <- function(
    stage,
    genes,
    view = 'sagittal',
    slices = NULL,
    annotation_level = 'custom_2',
    annotation_colors = many
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
    }

    possible_views <- c('sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice')
    if (!view %in% possible_views){
        stop(cat(
            paste0('"', view, '" is not a valid view argument. Try one of these:\n'),
            paste(possible_views, collapse = ', '), '\n'
        ))
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    voxel_meta <- DATA_LIST[[stage]]$row_meta

    if (view == 'slice'){
        slices <- seq(1, 40, 2)
        voxel_meta$voxel <- as.character(voxel_meta$voxel)
        expr_markers <- get_markers(voxel_mat, genes, scale=T) %>%
            mutate(voxel=as.character(voxel))
        expr_markers <- suppressMessages(right_join(expr_markers, voxel_meta))
        plot_df <- expr_markers %>%
            filter(x%in%slices) %>%
            mutate(slice=factor(x),
                   corr=expr,
                   group=factor(gene, levels=genes))
        p <- slice_plot(
            slice_df = plot_df,
            annotation_colors = annotation_colors,
            annotation_level = annotation_level,
            map_colors = gyrdpu
        )
        return(p)
    } else {
        if (view %in% c('coronal', 'x')){
            if (!is.null(slices)){
                voxel_meta <- filter(voxel_meta, x %in% slices)
            }
            f_plot <- function(x, g){
                ggplot(x, aes(z, y, fill=expr)) +
                    geom_tile() +
                    theme_bw() +
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
                ggplot(x, aes(x, z, fill=expr)) +
                    geom_tile() +
                    theme_bw() +
                    dr_theme +
                    labs(title=g) +
                    feature_fill_scale +
                    feature_theme
            }
        } else if (view %in% c('sagittal', 'z')){
            if (!is.null(slices)){
                voxel_meta <- filter(voxel_meta, z %in% slices)
            }
            f_plot <- function(x, g){
                ggplot(x, aes(x, y, fill=expr)) +
                    geom_tile() +
                    theme_bw() +
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
}


#' Plot gene expression across the mouse brain in 3D
#'
#' @import dplyr
#'
#' @export
#'
plot_expression_3d <- function(
    stage,
    gene,
    annotation_level = NULL,
    annotation_colors = many,
    expression_colors = inferno,
    sizes = c(10, 1000)
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    voxel_meta <- DATA_LIST[[stage]]$row_meta

    voxel_meta$voxel <- as.character(voxel_meta$voxel)
    expr_markers <- get_markers(voxel_mat, gene, scale=T) %>%
        mutate(voxel=as.character(voxel))
    expr_markers <- suppressMessages(right_join(expr_markers, voxel_meta)) %>%
        mutate(intensity=expr)
    plot_df <- expr_markers %>%
        mutate(z=-z+max(z)*2) %>%
        bind_rows(expr_markers, .id = 'hemisphere')
    p <- three_dim_plot(
        int_df = plot_df,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        intensity_colors = expression_colors,
        sizes = sizes
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
    ncol = NULL,
    nrow = NULL,
    sort = T,
    scale = T
){
    meta$voxel <- as.character(meta$voxel)
    expr_markers <- get_markers(expr_mat, markers, scale=scale) %>%
        mutate(voxel=as.character(voxel))
    expr_markers <- suppressMessages(right_join(expr_markers, meta))
    plots <- map(unique(markers), function(g){
        x <- expr_markers %>%
            filter(gene==g)
        if (sort){
            x <- arrange(x, expr)
        }
        plot(x, g)
    })
    p <- egg::ggarrange(plots = plots, ncol = ncol, nrow = nrow, draw = F)
    return(p)
}



#' 3D intensity plot
#'
three_dim_plot <- function(
    int_df,
    annotation_level = NULL,
    annotation_colors = many,
    intensity_colors = inferno,
    sizes = c(10, 1000)
){

    if (is.null(annotation_level)){
        p <- plotly::plot_ly(
            int_df,
            x=~x,
            y=~y,
            z=~z,
            size=~intensity,
            color=~intensity,
            sizes = sizes,
            colors=intensity_colors
        )
    } else {
        p <- plotly::plot_ly(
            int_df,
            x=~x,
            y=~y,
            z=~z,
            size=~intensity,
            color=int_df[[annotation_level]],
            sizes = sizes,
            colors=annotation_colors
        )
    }
    return(p)
}


