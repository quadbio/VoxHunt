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


#' Plot brain annotation
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#'
brain_annotation_plot <- function(stage = 'E13pt5',
                                  annotation_level = 'custom_2',
                                  annotation_colors = many,
                                  show_coordinates = F,
                                  show_legend = F,
                                  alpha = 0.5
){
    utils::data(voxel_meta, envir = environment())
    m <- filter(voxel_meta, age==stage)
    p <- ggplot(m, aes_string('x', 'y', fill=annotation_level)) +
        geom_tile(alpha=alpha)
    if (!show_coordinates){
        p <- p + theme_void()
    }
    if (!show_legend){
        p <- p + theme(legend.position = 'none')
    }
    return(p)
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

# single_feature_plot <- function(expr_mat, meta, gene,
#                                 plot=tsne_plot,
#                                 title=NULL,
#                                 sort=T,
#                                 scale=T){
#     meta$cell <- as.character(meta$cell)
#     expr_markers <- get_markers(expr_mat, gene, scale=scale) %>%
#         mutate(cell=as.character(cell)) %>%
#         right_join(meta)
#     x <- expr_markers %>%
#         filter(gene==gene)
#     if (sort){
#         x <- arrange(x, expr)
#     }
#     return(plot(x, gene))
# }
#
# heatmap_plot <- function(expr_mat, groups, genes,
#                          colors = blues_flat,
#                          group_colors = many,
#                          newpage = T
# ){
#
#     top_markers <- unique(genes)
#
#     top_marker_expr <- get_markers(expr_mat, top_markers) %>%
#         mutate(groups=factor(groups[match(.$cell, rownames(expr_mat))], levels=levels(factor(groups)))) %>%
#         group_by(groups) %>%
#         mutate(mean_expr=mean(expr)) %>%
#         ungroup() %>%
#         arrange(groups, mean_expr) %>%
#         mutate(cell=factor(cell, levels=unique(cell)),
#                gene=factor(gene, levels=top_markers))
#
#     p1 <- ggplot(top_marker_expr, aes(cell, gene, fill=expr)) +
#         scale_fill_gradientn(colors=colors) +
#         geom_raster() +
#         scale_y_discrete(expand=expand_scale(add = c(0,0))) +
#         theme(
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank(),
#             axis.ticks.x = element_blank(),
#             axis.ticks.y = element_blank(),
#             panel.border = element_blank(),
#             legend.position = 'none'
#         )
#     p2 <- ggplot(top_marker_expr, aes(cell, 1, fill=groups)) +
#         geom_raster() +
#         scale_fill_manual(values=group_colors) +
#         scale_y_discrete(expand=expand_scale(add = c(0,0)), labels=function(breaks) {rep_along(breaks, " ")}) +
#         theme(
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank(),
#             axis.ticks.x = element_blank(),
#             panel.border = element_blank(),
#             legend.position="bottom"
#         )
#     ggarrange(p1, p2, padding=0, heights = c(40, 2), newpage = newpage)
# }
