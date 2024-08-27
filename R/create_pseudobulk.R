#' Create Pseudo-bulked Assay data
#'
#' @param seurat_obj Seurat object
#' @param cell_type_var Chr name of metadata variable containing cell type annotation
#' @param grouping_vars Chr vector of metadata variables that define cell grouping strategy
#' @param split Lgl whether or not to split pseudo-bulk matrix by cell type
#'
#' @import Seurat
#' @import SeuratObject
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @return Data frame, or list of data frames if `split == TRUE`
#' @export
#'
#' @examples
#' create_pseudobulk(
#'     seurat_obj = Neurons,
#'     cell_type_var = "CellType",
#'     grouping_vars = c("Sex", "GEM"),
#'     split = TRUE
#' )
create_pseudobulk <- function(seurat_obj, cell_type_var, grouping_vars, split = FALSE) {
    cts <- AggregateExpression(
        object = seurat_obj,
        assays = "RNA",
        slot = "counts",
        group.by = c(cell_type_var, grouping_vars),
        return.seurat = FALSE
    )$RNA %>%
        t() %>%
        data.frame()

    if (split == TRUE) {
        # split by cell type
        splitRows <- gsub("_.*", "", rownames(cts))
        cts.split <- split.data.frame(cts,
            f = factor(splitRows)
        )

        # clip cell type from sample name, and return to gene x sample orientation
        cts.split.modified <- lapply(cts.split, function(x) {
            rownames(x) <- gsub("^([^_]*_){2}", "\\1", rownames(x))
            rownames(x) <- gsub("-", "_", rownames(x))
            x <- t(x) %>% data.frame()
            return(x)
        })

        return(cts.split.modified)
    }

    cts <- cts %>%
        t() %>%
        data.frame()
    return(cts)
}
