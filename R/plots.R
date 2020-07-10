#### ANNOTATION PLOTS ####

#' @import dplyr
#'
#' @param object String indicating the developmental stage from the ABA.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors A character vector with the color scale.
#' @param slices A numeric vector indicating the slices to plot.
#' @param show_coordinates Logical. Whether to show slice coordinates or not.
#' @param show_legend Logical. Whether to show a color legend or not.
#' @param alpha Float. The alpha value for grid tiles.
#'
#' @return An annotation plot.
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
    show_coordinates = FALSE,
    show_legend = FALSE,
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
    show_coordinates = FALSE,
    show_legend = FALSE,
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
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors A character vector with the color scale.
#' @param show_coordinates Logical. Whether to show slice coordinates or not.
#' @param show_legend Logical. Whether to show a color legend or not.
#' @param alpha Float. The alpha value for grid tiles.
#'
#' @return An annotation plot.
#'
annotation_plot <- function(
    annot_df,
    annotation_level = 'custom_2',
    annotation_colors = many,
    show_coordinates = FALSE,
    show_legend = FALSE,
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

#### PLOTS FOR VOXEL MAPS ####

#' @import ggplot2
#' @import dplyr
#'
#' @param view String indicating the perspective to show. Valid values are
#' 'sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'.
#' @param slices A numeric vector indicating the slices to plot.
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param map_colors Color map for correlation values.
#' @param slices A numeric vector indicating the slices to plot.
#' @param show_coordinates Logical. Whether to show slice coordinates or not.
#' @param show_legend Logical. Whether to show a color legend or not.
#' @param return_plot Logical. Return plot.
#' @param ... Other arguments passed to egg::ggarrange().
#'
#' @return A similarity map.
#'
#' @rdname plot_map
#' @export
#' @method plot_map VoxelMap
#'
plot_map.VoxelMap <- function(
    object,
    view = c('sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'),
    slices = NULL,
    groups = NULL,
    annotation_level = 'custom_2',
    annotation_colors = many,
    map_colors = gyrdpu_flat,
    show_coordinates = FALSE,
    show_legend = FALSE,
    return_plot = FALSE,
    ...
){
    view <- match.arg(view)

    plot_df <- summarize_groups(object, groups)

    if (view == 'slice'){
        if (is.null(slices)){
            slices <- seq(1, 40, 2)
        }
        plot_df <- plot_df %>%
            filter(x%in%slices) %>%
            mutate(slice=factor(x))

        p <- slice_plot(
            slice_df = plot_df,
            annotation_colors = annotation_colors,
            annotation_level = annotation_level,
            map_colors = map_colors,
            ...
        )
    } else {
        p <- mapping_plot(
            map_df = plot_df,
            slices = slices,
            view = view,
            map_colors = map_colors,
            ...
        )
    }
    return(p)
}


#' @import ggplot2
#' @import dplyr
#'
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param group_show The nameof the group to show in the plot.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param map_colors Color map for correlation values.
#' @param sizes The size range for points.
#' @param both_hemispheres Logical. Whether to plot both hemispheres (TRUE) of only one (FALSE).
#' @param ... Other arguments passed to plotly::plot_ly().
#'
#' @return A similarity map in 3D.
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
    map_colors = inferno_flat,
    sizes = c(1, 1000),
    both_hemispheres = TRUE,
    ...
){

    plot_df <- summarize_groups(object, groups)

    if (is.null(show_group) & is.null(groups)){
        show_group <- levels(factor(object$cell_meta$group))[1]
    } else if (is.null(show_group) & !is.null(groups)){
        show_group <- levels(factor(groups))[1]
    }

    plot_df <- filter(plot_df, group==show_group) %>%
        mutate(intensity=corr)
    if (both_hemispheres){
        plot_df <- plot_df %>%
            mutate(z=-z+max(z)*2) %>%
            bind_rows(plot_df, .id = 'hemisphere')
    }
    p <- three_dim_plot(
        int_df = plot_df,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        intensity_colors = map_colors,
        sizes = sizes,
        ...
    )
    return(invisible())
}


#' Plot correlation to brain slices and structure annotation
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param map_colors Color map for correlation values.
#' @param ... Other arguments passed to egg::ggarrange().
#'
#' @return A similarity map with one row per slice and one
#' column per cell group.
#'
slice_plot <- function(
    slice_df,
    annotation_level = 'custom_2',
    annotation_colors = many,
    map_colors = gyrdpu,
    ...
){
    annot <- ggplot(slice_df, aes_string('z', 'y', fill = annotation_level)) +
        geom_tile() +
        theme_void() +
        scale_fill_manual(values = annotation_colors) +
        facet_grid(slice~., switch = 'both') +
        theme(
            legend.position = 'right',
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
                strip.text = element_blank(),
                legend.position = 'none',
                plot.title = element_text(angle = 60, hjust = 0, vjust = 0, size = 10)
                ) +
            labs(title = g)
        return(p)
    })
    plots <- c(list(annot), plots)
    out_plot <- patchwork::wrap_plots(plots, nrow = 1, ...) +
        patchwork::plot_layout(guides = 'collect')
    return(out_plot)
}


#' Plot correlation to entire brain or selected slices.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param view String indicating the perspective to show. Valid values are
#' 'sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'.
#' @param slices A numeric vector indicating the slices to plot.
#' @param map_colors Color map for correlation values.
#' @param ... Other arguments passed to egg::ggarrange().
#'
#' @return A similarity map.
#'
mapping_plot <- function(
    map_df,
    view = 'sagittal',
    slices = NULL,
    map_colors = rdpu,
    ...
){
    plots <- map(levels(factor(map_df$group)), function(g){
        plot_df <- filter(map_df, group==g)
        if (view %in% c('coronal', 'x')){
            if (!is.null(slices)){
                plot_df <- filter(plot_df, x %in% slices)
            }
            coords <- c('y', 'z')
        } else if (view %in% c('traverse', 'y')){
            if (!is.null(slices)){
                plot_df <- filter(plot_df, y %in% slices)
            }
            coords <- c('x', 'z')
        } else {
            if (!is.null(slices)){
                plot_df <- filter(plot_df, z %in% slices)
            }
            coords <- c('x', 'y')
        }
        p <- ggplot(plot_df, aes_string(coords[1], coords[2], fill = 'corr', alpha = 'corr')) +
            geom_tile() +
            theme_bw() +
            theme_void() +
            theme(
                plot.title = element_text(size = 8)
                ) +
            scale_fill_gradientn(colors = map_colors) +
            scale_alpha_continuous(guide = FALSE) +
            labs(title = g)
        return(p)
    })
    return(patchwork::wrap_plots(plots, ...))
}


#### STRUCTURE SIMILARITY PLOTS ####

#' @import ggplot2
#' @import dplyr
#'
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param cluster Logical. Perform hierarchical clustering on rows and columns.
#' @param return_plot Logical. Return plot.
#'
#' @return A plot showing average correlation to each brain structure.
#'
#' @rdname plot_structure_similarity
#' @export
#' @method plot_structure_similarity VoxelMap
#'
plot_structure_similarity.VoxelMap <- function(
    object,
    groups = NULL,
    annotation_level = 'custom_3',
    annotation_colors = blues_flat,
    cluster = TRUE,
    return_plot = FALSE,
    include_unannotated = FALSE,
    ...
){
    annot_df <- object$voxel_meta %>%
        distinct_('custom_2', annotation_level) %>%
        rename('struct'=annotation_level, 'annot'='custom_2')
    if (annotation_level == 'custom_2'){
        annot_df$struct <- annot_df$annot
    }
    annot_df <- suppressMessages(inner_join(annot_df, struct_custom2)) %>%
        mutate(symbol=factor(symbol, levels=struct_symbol_custom2)) %>%
        filter(!is.na(struct))

    plot_df <- summarize_groups(object, groups) %>%
        group_by_('group', annotation_level) %>%
        summarize(corr=mean(corr)) %>%
        ungroup() %>%
        group_by(group) %>%
        mutate(corr=zscale(corr)) %>%
        ungroup() %>%
        rename('struct'=annotation_level) %>%
        mutate(group=as.character(group), struct=as.character(struct))

    if (cluster){
        clust_df <- spread(plot_df, struct, corr) %>%
            as_matrix()
        group_clust <- as.dendrogram(hclust(dist(clust_df), 'ward.D2'))
        struct_clust <- as.dendrogram(hclust(dist(t(clust_df)), 'ward.D2'))
        plot_df <- suppressMessages(left_join(plot_df, annot_df)) %>%
            mutate(
                group = factor(group, levels=labels(group_clust)),
                struct = factor(struct, levels=labels(struct_clust))
            )
    } else {
        plot_df <- suppressMessages(left_join(plot_df, annot_df)) %>%
            arrange(symbol) %>%
            mutate(struct = factor(struct, levels=unique(struct)))
    }

    if (!include_unannotated){
        plot_df <- filter(plot_df, !stringr::str_detect(annot, 'telencephalic|transition'))
    }

    p1 <- ggplot(plot_df, aes(struct, group, fill=corr)) +
        geom_tile() +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        scale_fill_gradientn(colors=blues_flat) +
        theme_article() +
        theme(
            axis.text.x = element_text(angle=45, hjust=1)
        ) +
        labs(x = 'Structure', y = 'Group', fill = 'Scaled\ncorrelation')

    p2 <- ggplot(plot_df, aes(struct, fill=symbol)) +
        geom_bar(position='fill', width=1) +
        theme_void() +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        scale_fill_manual(values=struct_colors_custom2) +
        theme(
            legend.position = 'top'
        ) +
        labs(fill='Structure\nannotation')

    p <- egg::ggarrange(p2, p1, heights=c(1,10), ncol = 1, ...)
    if (return_plot){
        return(p)
    } else {
        return(invisible())
    }
}


#### EXPRESSION PLOTS ####

#' Plot gene expression across the mouse brain
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param stage The developmental stage to plot.
#' @param genes A character vector with the genes to plot.
#' @param view String indicating the perspective to show. Valid values are
#' 'sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'.
#' @param slices A numeric vector indicating the slices to plot.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param return_plot Logical. Return plot.
#' @param ... Other arguments passed to egg::ggarrange().
#'
#' @return A gene expression plot.
#'
#' @export
#'
plot_expression <- function(
    stage,
    genes,
    view = 'sagittal',
    slices = NULL,
    annotation_level = 'custom_2',
    annotation_colors = many,
    return_plot = FALSE,
    ...
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
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
        if (is.null(slices)){
            slices <- seq(1, 40, 2)
        }
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
            map_colors = gyrdpu,
            ...
        )
    } else {
        if (view %in% c('coronal', 'x')){
            if (!is.null(slices)){
                voxel_meta <- filter(voxel_meta, x %in% slices)
            }
            coords <- c('z', 'y')
        } else if (view %in% c('traverse', 'y')) {
            if (!is.null(slices)){
                voxel_meta <- filter(voxel_meta, y %in% slices)
            }
            coords <- c('x', 'z')
        } else if (view %in% c('sagittal', 'z')){
            if (!is.null(slices)){
                voxel_meta <- filter(voxel_meta, z %in% slices)
            }
            coords <- c('x', 'y')
        }
        f_plot <- function(x, g){
            ggplot(x, aes_string(coords[1], coords[2], fill='expr')) +
                geom_tile() +
                theme_bw() +
                theme_void() +
                dr_theme +
                labs(title=g) +
                feature_fill_scale +
                feature_theme
        }
        p <- feature_plot(
            expr_mat = voxel_mat,
            meta = voxel_meta,
            markers = genes,
            plot = f_plot,
            ...
        )
    }
    if (return_plot){
        return(p)
    } else {
        return(invisible())
    }
}


#' Plot gene expression across the mouse brain in 3D
#'
#' @import dplyr
#'
#' @param stage The developmental stage to plot.
#' @param genes A character vector with the genes to plot.
#' @param view String indicating the perspective to show. Valid values are
#' 'sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'.
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param expression_colors Colors for expression scale.
#' @param sizes The size range for points.
#' @param both_hemispheres Logical. Whether to plot both hemispheres (TRUE) of only one (FALSE).
#' @param ... Other arguments passed to plotly::plot_ly().
#'
#' @return A gene expression plot in 3D.
#'
#' @export
#'
plot_expression_3d <- function(
    stage,
    gene,
    annotation_level = NULL,
    annotation_colors = many,
    expression_colors = inferno,
    sizes = c(10, 1000),
    both_hemispheres = TRUE,
    ...
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    voxel_meta <- DATA_LIST[[stage]]$row_meta

    voxel_meta$voxel <- as.character(voxel_meta$voxel)
    expr_markers <- get_markers(voxel_mat, gene, scale=TRUE) %>%
        mutate(voxel=as.character(voxel))
    plot_df <- suppressMessages(right_join(expr_markers, voxel_meta)) %>%
        mutate(intensity=expr)
    if (both_hemispheres){
        plot_df <- plot_df %>%
            mutate(z=-z+max(z)*2) %>%
            bind_rows(plot_df, .id = 'hemisphere')
    }
    p <- three_dim_plot(
        int_df = plot_df,
        annotation_level = annotation_level,
        annotation_colors = annotation_colors,
        intensity_colors = expression_colors,
        sizes = sizes,
        ...
    )
    return(invisible())
}


#' Feature plot
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param expr_mat A sample x gene expression matrix.
#' @param meta Metadata with embedding coordinates.
#' @param markers A character vector with genes to display.
#' @param plot Function to use for generating individual plots.
#' @param sort Logical. Whether to sort highest expressing cells up.
#' @param scale Logical. Whether to scale the expression values.
#' @param ... Other arguments passed to egg::ggarrange().
#'
#' @return A panel of feature plots for each gene.
#'
feature_plot <- function(
    expr_mat,
    meta,
    markers,
    plot = xy_plot,
    sort = TRUE,
    scale = TRUE,
    ...
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
    return(egg::ggarrange(plots = plots, ...))
}

#### 3D PLOT ####

#' General 3D plotting function
#'
#' @param annotation_level The structure annotation level to color code.
#' @param annotation_colors Color map for structure annotation.
#' @param intnesity_colors Colors for intensity scale.
#' @param sizes The size range for points.
#' @param both_hemispheres Logical. Whether to plot both hemispheres (TRUE) of only one (FALSE).
#' @param ... Other arguments passed to plotly::plot_ly().
#'
#' @return A three-dimensional plot.
#'
three_dim_plot <- function(
    int_df,
    annotation_level = NULL,
    annotation_colors = many,
    intensity_colors = inferno,
    sizes = c(10, 1000),
    ...
){
    if (is.null(annotation_level)){
        p <- plotly::plot_ly(
            int_df,
            x = ~x,
            y = ~y,
            z = ~z,
            size = ~intensity,
            color = ~intensity,
            sizes = sizes,
            colors = intensity_colors,
            type = 'scatter3d',
            mode = 'markers',
            ...
        )
    } else {
        p <- plotly::plot_ly(
            int_df,
            x = ~x,
            y = ~y,
            z = ~z,
            size = ~intensity,
            color = int_df[[annotation_level]],
            sizes = sizes,
            colors = annotation_colors,
            type = 'scatter3d',
            mode = 'markers',
            ...
        )
    }
    axis <- list(
        showgrid = F
    )
    scene <- list(
        aspectmode='data',
        xaxis = axis,
        yaxis = axis,
        zaxix = axis
    )
    p <- plotly::layout(p, scene=scene)
    suppressWarnings(print(p))
}


#### PLOTS FOR BRAINSPAN MAPS ####
#' @import ggplot2
#' @import dplyr
#'
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param annotation_level The structure annotation level to summarize to.
#'
#' @return A similarity map.
#'
#' @rdname plot_map
#' @export
#' @method plot_map ReferenceMap
#'
plot_map.ReferenceMap <- function(
    object,
    groups = NULL,
    annotation_level = c('structure_group', 'structure_name', 'structure_acronym'),
    map_colors = blues_flat
){
    annotation_level <- match.arg(annotation_level)
    if (annotation_level == 'structure_name'){
        annotation_level <- 'structure_acronym'
    }

    plot_df <- summarize_groups(object, groups) %>%
        group_by_('group', annotation_level) %>%
        dplyr::summarize(corr=mean(corr)) %>%
        ungroup() %>%
        group_by(group) %>%
        mutate(corr=zscale(corr))

    p <- ggplot(plot_df, aes_string('group', annotation_level, fill='corr')) +
        geom_tile() +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        scale_fill_gradientn(colors=blues_flat) +
        theme_article() +
        theme(
            axis.text.x = element_text(angle=45, hjust=1)
        ) +
        labs(x = 'Group', y = 'Structure', fill = 'Scaled\ncorrelation')
    return(p)
}

#### UTILS ####
#' @import ggplot2
#'
#' @rdname no_legend
#' @export
#'
no_legend <- function(){
    theme(
        legend.position = 'none'
    )
}
