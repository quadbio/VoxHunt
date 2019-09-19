#' Aggregate matrix over groups
#'
#' @import Matrix
#' @import dplyr
#'
aggregate_matrix <- function(x,
                             groups = NULL,
                             fun = colMeans
                             ){
    if (!is.null(groups)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            fun(x[which(groups==g), ])
        })
        rownames(agg_mat) <- colnames(x)
    } else {
        agg_mat <- as.matrix(fun(x))
        colnames(agg_mat) <- ' '
    }
    agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
    return(agg_mat)
}



# extract_legend <- function(x){
#     tmp <- ggplot::ggplot_gtable(ggplot::ggplot_build(x))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     return(legend)}
