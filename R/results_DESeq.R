#' Get DESeq Results
#'
#' @param DESeq_obj DESeq object
#' @param seurat_obj Seurat object
#' @param Contrast Chr vector of contrast in the format c("variable", "level1", "level2")
#'
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import DESeq2
#'
#' @return Data frame of DESeq results
#' @export
#'
#' @examples
#' # results_DESeq(
#' #   seurat_obj = Neurons,
#' #   DESeq_obj = DESeq_obj,
#' #   Contrast = c("Treatment", "Mor", "Sal")
#' # )
#'
results_DESeq <- function(DESeq_obj, seurat_obj, Contrast, cell_type_var, cell_type_level) {
  # Extract results
  res <- results(DESeq_obj, contrast = Contrast) %>%
    data.frame() %>%
    rownames_to_column(var = "GeneName")
  # Extract counts table
  counts_mat <- seurat_obj@assays[["RNA"]]@counts
  # Replace any hyphens in gene names with ".", to match results table
  rownames(counts_mat) <- gsub("-", ".", rownames(counts_mat))
  # DESeq2 will prepend "X" to genes with names that start with numbers
  # consequently, this mucks up a bunch of stuff so we'll get ahead of it here
  # Identify genes that start with a number
  numeric_genes <- grep("^[0-9]", rownames(counts_mat), value = TRUE)
  # modify the corresponding gene names in DESeq2 results
  modified_gene_names <- paste0("X", numeric_genes)
  # create a mapping between the modified and original gene names
  gene_mapping <- setNames(numeric_genes, modified_gene_names)
  res$GeneName <- ifelse(res$GeneName %in% names(gene_mapping),
                         gene_mapping[res$GeneName],
                         res$GeneName
  )
  # Format results table
  cells.1 <- seurat_obj@meta.data %>%
    filter(!!sym(cell_type_var) == cell_type_level) %>%
    filter(!!sym(Contrast[1]) == Contrast[2]) %>%
    row.names()
  cells.2 <- setdiff(x = seurat_obj@meta.data %>%
                       filter(!!sym(cell_type_var) == cell_type_level) %>%
                       row.names(),
                     y = cells.1)
  res$Pct_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.1, drop = FALSE] > 0) / length(cells.1)) * 100
  res$Pct_Other_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.2, drop = FALSE] > 0) / length(cells.2)) * 100
  res$CellType <- cell_type_level
  # order column names
  res <- res %>%
    select(CellType, GeneName, baseMean, log2FoldChange, lfcSE, stat, pvalue,
           padj, Pct_Expressing, Pct_Other_Expressing
    )
  return(res)
}
