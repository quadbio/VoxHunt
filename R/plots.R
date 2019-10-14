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
#' @import ggplot2
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
    show_coordinates = F,
    show_legend = F,
    alpha = 0.5
){
    utils::data(voxel_meta, envir = environment())
    age <- stage_name(object)
    m <- filter(voxel_meta, stage == age)
    p <- ggplot(m, aes_string('x', 'y', fill = annotation_level)) +
        geom_tile(alpha = alpha) +
        scale_x_continuous(breaks = seq(0, 100, 2)) +
        scale_y_continuous(breaks = seq(-100, 100, 2)) +
        labs(title = age)
    if (!show_coordinates){
        p <- p + theme_void()
    }
    if (!show_legend){
        p <- p + theme(legend.position = 'none')
    }
    return(p)
}


#' Plot brain structure annotation
#'
#' @import ggplot2
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
    show_coordinates = F,
    show_legend = F,
    alpha = 0.5
){
    m <- object$voxel_meta
    age <- unique(m$stage)
    p <- ggplot(m, aes_string('x', 'y', fill = annotation_level)) +
        geom_tile(alpha = alpha) +
        scale_x_continuous(breaks = seq(0, 100, 2)) +
        scale_y_continuous(breaks = seq(-100, 100, 2)) +
        labs(title = age)
    if (!show_coordinates){
        p <- p + theme_void()
    }
    if (!show_legend){
        p <- p + theme(legend.position = 'none')
    }
    return(p)
}


#' Plot correlation to brain slices and structure annotation
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
slice_plot <- function(slice_df,
                       struct_colors = many,
                       corr_colors = gyrdpu,
                       newpage = T){

    annot <- ggplot(slice_df, aes(z, y, fill = struct)) +
        geom_tile() +
        theme_void() +
        scale_fill_manual(values = struct_colors) +
        facet_grid(slice~.) +
        theme(
            legend.position = 'none',
            strip.text = element_blank()
        )

    plots <- purrr::map(levels(factor(slice_df$cluster)), function(x){
        plot_df <- slice_df %>%
            {.[.$cluster==x, ]} %>%
            dplyr::group_by(z, y, slice) %>%
            dplyr::summarize(corr = mean(corr, na.rm = T))

        p <- ggplot(plot_df, aes(-z, y, fill = corr, alpha = corr)) +
            geom_tile() +
            theme_void() +
            scale_fill_gradientn(colors = corr_colors) +
            scale_alpha_continuous(range = c(0.5, 1)) +
            facet_grid(slice ~ .) +
            theme(
                legend.position = 'none',
                strip.text = element_blank(),
                plot.title = element_text(angle = 45, hjust = 0, vjust = 0, size = 8)
                ) +
            labs(title = x)
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
#' @export
#'
brain_plot <- function(corr_df,
                       corr_colors = rdpu,
                       newpage = T
){
    plots <- map(levels(factor(corr_df$cluster)), function(x){
        plot_df <- corr_df[corr_df$cluster==x, ]
        p <- ggplot(plot_df, aes(x, y, fill = corr)) +
            geom_tile() +
            theme_void() +
            theme(
                legend.position = 'none',
                plot.title = element_text(angle = 45, hjust = 0, vjust = 0, size = 8)
                ) +
            scale_fill_gradientn(colors = corr_colors) +
            labs(title = x)
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
#'
brain_expression_plot <- function(expr_mat,
                                  stage = 'E13pt5',
                                  genes = NULL

){
    utils::data(voxel_meta, envir = environment())
    m <- filter(voxel_meta, age==stage)
    feature_plot(expr_mat = expr_mat, meta = m, markers = genes)
}

#' Feature plot
#'
#' @import ggplot2
#' @import dplyr
#'
feature_plot <- function(expr_mat, meta, markers,
                         plot = xy_plot,
                         title = NULL,
                         ncol = NULL,
                         nrow = NULL,
                         sort = T,
                         scale = T){
    meta$cell <- as.character(meta$voxel)
    expr_markers <- get_markers(expr_mat, markers, scale=scale) %>%
        mutate(cell=as.character(cell)) %>%
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

#' Plotting xy coordinates
#'
#' @import ggplot2
#' @import dplyr
#'
xy_plot <- function(x, g){
    ggplot(x, aes(x, y, fill=expr, alpha=expr)) +
        geom_tile(alpha=point_alpha) +
        dr_theme +
        labs(title=g) +
        feature_fill_scale +
        feature_theme
}


