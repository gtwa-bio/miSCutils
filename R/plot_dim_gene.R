#' Plot by cell gene expression
#'
#' Plot gene expression (in log counts) for a particular gene within cells in
#' a reduced dimensionality plot. This function is borrowed in large part from
#' the iSEE reduced dimensionality plot panel.
#'
#' @param SCE Single Cell Experiment object
#' @param Dim String. Name of dimensionality reduciton to plot
#' @param Gene String. Name of gene to plot
#' @param Exp_Order Boolean. If TRUE, cells will be plotted in order of count values
#'
#' @import SingleCellExperiment
#' @importFrom iSEE ExperimentColorMap synchronizeAssays
#' @import ggplot2
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # marker_plot(SCE = se,
#' #             Dim = "PC12UMAP",
#' #             Gene = "Th",
#' #             Exp_Order = TRUE
#' # )
plot_dim_geme <- function(SCE, Dim, Gene, Exp_Order = FALSE){
  # plot setup
  colormap <- iSEE::ExperimentColorMap()
  colormap <- iSEE::synchronizeAssays(colormap, SCE)

  # select reduction
  red.dim <- reducedDim(SCE, Dim);
  plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(SCE))

  # select assay and feature to plot
  plot.data$ColorBy <- assay(SCE, "logcounts")[Gene, ]

  if(Exp_Order == TRUE){
    # Order plot.data by ColorBy
    plot.data <- plot.data[order(plot.data$ColorBy), ]
  } else {
    plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE]
  }

  # plot
  dot.plot <- ggplot() +
    geom_point(aes(x=X, y=Y, color=ColorBy), alpha=1, plot.data, size=1) +
    labs(x="UMAP 1", y="UMAP 2", color=paste0(Gene, "\n(logcounts)"), title=Gene) +
    coord_cartesian(xlim=range(plot.data$X, na.rm=TRUE),
                    ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_color_gradientn(colors=c("lightgrey", "blue"), na.value='grey50', limits=range(plot.data$ColorBy, na.rm=TRUE)) +
    theme_bw() +
    theme(legend.position='bottom', legend.box='vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
          axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))

  # return plot
  dot.plot
}
