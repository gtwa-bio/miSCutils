#' Plot cluster integration as stacked bar plot
#'
#' @param seurat_obj
#' @param cluster_col
#' @param treatment_col
#' @param gem_col
#' @param sex_col
#'
#' @import Seurat
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#'
#' @return
#' @export
#'
#' @examples
plot_integration_bar <- function(seurat_obj, cluster_col, treatment_col, gem_col, sex_col) {
  # Treatment Group
  T_count_table <- table(seurat_obj@meta.data[[cluster_col]], seurat_obj@meta.data[[treatment_col]])
  T_count_mtx   <- as.data.frame.matrix(T_count_table)
  T_count_mtx$cluster <- rownames(T_count_mtx)
  T_melt_mtx    <- melt(T_count_mtx)
  T_melt_mtx$cluster <- factor(T_melt_mtx$cluster, levels = levels(seurat_obj@meta.data[[cluster_col]]))

  # will be the same for all integration features, will only calculate once, begin
  cluster_size   <- aggregate(value ~ cluster, data = T_melt_mtx, FUN = sum)
  cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(seurat_obj@meta.data[[cluster_col]]))
  # end

  colnames(T_melt_mtx)[2] <- "Treatment"

  # GEM well
  G_count_table <- table(seurat_obj@meta.data[[cluster_col]], seurat_obj@meta.data[[gem_col]])
  G_count_mtx   <- as.data.frame.matrix(G_count_table)
  G_count_mtx$cluster <- rownames(G_count_mtx)
  G_melt_mtx    <- melt(G_count_mtx)
  G_melt_mtx$cluster <- factor(G_melt_mtx$cluster)
  G_melt_mtx$cluster <- factor(G_melt_mtx$cluster,levels = levels(seurat_obj@meta.data[[cluster_col]]))
  colnames(G_melt_mtx)[2] <- "GEM"

  # Predicted Sex
  S_count_table <- table(seurat_obj@meta.data[[cluster_col]], seurat_obj@meta.data[[sex_col]])
  S_count_mtx   <- as.data.frame.matrix(S_count_table)
  S_count_mtx$cluster <- rownames(S_count_mtx)
  S_melt_mtx    <- melt(S_count_mtx)
  S_melt_mtx$cluster <- factor(S_melt_mtx$cluster)
  S_melt_mtx$cluster <- factor(S_melt_mtx$cluster,levels = levels(seurat_obj@meta.data[[cluster_col]]))
  colnames(S_melt_mtx)[2] <- "Predicted.Sex"

  p1 <- ggplot(T_melt_mtx,aes(x=cluster,y=value,fill=Treatment)) +
    geom_bar(position="fill", stat="identity") +
    theme_bw() +
    coord_flip() +
    labs(y = "", x = "Cluster", title = "Treatment Group") +
    theme(legend.position="bottom", legend.direction = "vertical")

  p2 <- ggplot(G_melt_mtx,aes(x=cluster,y=value,fill=GEM)) +
    geom_bar(position="fill", stat="identity") +
    theme_bw() +
    coord_flip() +
    labs(y = "", x = "", title = "GEM Well") +
    theme(legend.position="bottom",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank())

  p3 <- ggplot(S_melt_mtx,aes(x=cluster,y=value,fill=Predicted.Sex)) +
    geom_bar(position="fill", stat="identity") +
    theme_bw() +
    coord_flip() +
    labs(y = "", x = "", title = "Predicted Sex") +
    theme(legend.position="bottom", legend.direction = "vertical",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank())

  p4 <- ggplot(cluster_size, aes(x=cluster, y=value)) +
    geom_bar(position="dodge", stat="identity",fill = "grey60") +
    theme_bw() +
    scale_y_log10() +
    coord_flip() +
    labs(y = "", x = "", title = "Cell Count") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank())

  combined_plot <- p1 + p2 + p3 + p4 + plot_layout(widths = c(5,5,5,2.5))

  return(combined_plot)
}
