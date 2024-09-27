#' Plot cluster integration by LISI as ridge plot
#'
#' @param seurat_obj
#' @param cluster_col
#' @param treatment_col
#' @param gem_col
#' @param sex_col
#' @param reduction
#'
#' @import Seurat
#' @import lisi
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#'
#' @return
#' @export
#'
#' @examples
plot_integration_lisi <- function(seurat_obj, cluster_col, treatment_col, gem_col, sex_col, reduction = "umap") {
    ## take an integrated Seurat object, plot distributions of LISI index for treatment group, GEM well, and predicted sex

    # Calculate for groups
    ## Treatment group
    seurat_obj@meta.data$Treatment.LISI <- compute_lisi(
        X = seurat_obj@reductions[[reduction]]@cell.embeddings,
        meta_data = seurat_obj@meta.data,
        label_colnames = treatment_col
    )[[treatment_col]]

    ## GEM well
    seurat_obj@meta.data$GEM.LISI <- compute_lisi(
        X = seurat_obj@reductions[[reduction]]@cell.embeddings,
        meta_data = seurat_obj@meta.data,
        label_colnames = gem_col
    )[[gem_col]]

    ## predicted sex
    seurat_obj@meta.data$Sex.LISI <- compute_lisi(
        X = seurat_obj@reductions[[reduction]]@cell.embeddings,
        meta_data = seurat_obj@meta.data,
        label_colnames = sex_col
    )[[sex_col]]

    # Plot ridges
    Treatment_lisi_ridge <- ggplot(
        seurat_obj@meta.data,
        aes(
            x = Treatment.LISI,
            y = !!sym(cluster_col),
            fill = !!sym(cluster_col)
        )
    ) +
        theme_bw() +
        ggridges::geom_density_ridges() +
        labs(x = "LISI", y = "Cluster", title = "Treatment Group") +
        NoLegend()

    GEM_lisi_ridge <- ggplot(
        seurat_obj@meta.data,
        aes(
            x = GEM.LISI,
            y = !!sym(cluster_col),
            fill = !!sym(cluster_col)
        )
    ) +
        ggridges::geom_density_ridges() +
        theme_bw() +
        labs(x = "LISI", y = "", title = "GEM Well") +
        NoLegend() +
        theme(
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
        )

    Sex_lisi_ridge <- ggplot(
        seurat_obj@meta.data,
        aes(
            x = Sex.LISI,
            y = !!sym(cluster_col),
            fill = !!sym(cluster_col)
        )
    ) +
        ggridges::geom_density_ridges() +
        theme_bw() +
        labs(x = "LISI", y = "", title = "Predicted Sex") +
        NoLegend() +
        theme(
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
        )

    combined.plot <- Treatment_lisi_ridge + GEM_lisi_ridge + Sex_lisi_ridge + plot_layout(widths = c(3, 3, 3), heights = c(5, 5, 5), ncol = 3, nrow = 1)
    return(combined.plot)
}
