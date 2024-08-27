#' Conduct Pseudo-bulked Differential Expression Within Cell Types
#'
#' This function conducts our standard workflow for pseudo-bulked differential
#' gene expression testing within cell types using DESeq2. For a more modular
#' and interactive workflow, see the individual `create_pseudobulk()`,
#' `create_metadata()`, `create_DESeq()`, `DESeq()`, and `results_DESeq()` functions.
#' If running this analysis for the first time, we'd recommend running each step
#' interactively so that you can inspect the model fit quality.
#'
#' @param seurat_obj Seurat object
#' @param cell_type_var Chr name of metadata variable containing cell type annotation
#' @param grouping_vars Chr vector of metadata variables that define cell grouping strategy
#' @param counts_filter Int minimum count to retain gene
#' @param Design Formula for model to fit
#' @param Contrast Chr vector of contrast in the format c("variable", "level2", "level1")
#' @param ... Additional arguments to pass to DESeq2
#'
#' @import Seurat
#' @import SeuratObject
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import DESeq2
#'
#' @return Data frame of DESeq results
#' @export
#'
#' @examples
#' within_celltype_DEG(
#'     seurat_obj = Neurons,
#'     cell_type_var = "CellType",
#'     grouping_vars = c("Sex", "GEM"),
#'     Design = ~Treatment,
#'     Contrast = c("Treatment", "Mor", "Sal")
#' )
within_celltype_DEG <- function(seurat_obj,
    cell_type_var,
    grouping_vars,
    counts_filter = 5,
    Design,
    Contrast,
    ...) {
    # Create pseudobulk count data
    cts <- create_pseudobulk(
        seurat_obj = seurat_obj,
        cell_type_var = cell_type_var,
        grouping_vars = grouping_vars,
        split = TRUE
    )

    # Create pseudobulk metadata
    metadata <- create_metadata(
        seurat_obj = seurat_obj,
        sample_id_vars = grouping_vars,
        metadata_vars = all.vars(Design),
        cell_type_var = cell_type_var,
        split = TRUE
    )

    # Create DESeq object
    DESeq_list <- create_DESeq(
        pseudobulk_counts = cts,
        pseudobulk_metadata = metadata,
        counts_filter = counts_filter,
        Design = Design,
        split = TRUE
    )

    # Fit DESeq object with model
    DESeq_list <- map(
        DESeq_list,
        DESeq,
        ...
    )

    # Extract DEG results
    res <- map(DESeq_list, results_DESeq,
        seurat_obj = seurat_obj,
        Contrast = Contrast
    )

    # Return results
    return(res)
}
