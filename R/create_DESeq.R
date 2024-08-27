#' Create DESeq object or list of DESeq objects
#'
#' @param pseudobulk_counts Data frame of gene by sample counts
#' @param pseudobulk_metadata Data frame of sample metadata
#' @param counts_filter Int minimum count to retain gene
#' @param Design Formula for model to fit
#' @param split Lgl is pseudobulk data split by celltype
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import DESeq2
#'
#'
#' @return DESeq object, or list of DESeq objects if `split == TRUE`
#' @export
#'
#' @examples
#' # create_DESeq(
#' #     pseudobulk_counts = counts.list,
#' #   pseudobulk_metadata = metadata.list,
#' #   counts_filter = 5,
#' #   Design = ~Treatment,
#' #   split = TRUE
#' # )
#'
create_DESeq <- function(pseudobulk_counts, pseudobulk_metadata, counts_filter, Design, split = FALSE) {
    if (split == TRUE) {
        # Make DESeq object for each cell type
        dds <- lapply(names(pseudobulk_counts), function(celltype) {
            celltype.dds <- DESeqDataSetFromMatrix(
                countData = pseudobulk_counts[[celltype]],
                colData = pseudobulk_metadata[[celltype]],
                design = Design
            )
            # Filter low count genes
            keep <- rowSums(counts(celltype.dds)) >= counts_filter
            celltype.dds <- celltype.dds[keep, ]
            # return cell type dds object
            return(celltype.dds)
        })
    } else {
        # Make DESeq object
        dds <- DESeqDataSetFromMatrix(
            countData = pseudobulk_counts,
            colData = pseudobulk_metadata,
            design = Design
        )
        # Filter low count genes
        keep <- rowSums(counts(dds)) >= counts_filter
        dds <- dds[keep, ]
    }
    # Return final object
    return(dds)
}
