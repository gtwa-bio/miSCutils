#' Create Pseudo-bulked metadata table
#'
#' @param seurat_obj Seurat object
#' @param sample_id_vars Chr vector of metadata fields to construct sample ID from
#' @param metadata_vars Chr vector of metadata fields to include in metadata output
#' @param cell_type_var Chr name of metadata field including cell type annotation
#' @param split Lgl weather or not metadata table should be split by cell type
#'
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @return Data frame, list of data frames if `split == TRUE`
#' @export
#'
#' @examples
#' # create_metadata(
#' #    seurat_obj = Neurons,
#' #   sample_id_vars = c("Sex", "GEM"),
#' #   metadata_vars = c("Sex", "GEM", "Treatment"),
#' #   cell_type_var = "CellType",
#' #   split = TRUE
#' # )
#'
create_metadata <- function(seurat_obj, sample_id_vars, metadata_vars, cell_type_var, split = FALSE) {
    # Extract metadata from the Seurat object
    metadata <- seurat_obj@meta.data

    # Create sample ID by combining specified columns
    if (split == TRUE) {
        metadata$sample.id <- apply(metadata[, sample_id_vars], 1, paste, collapse = "_")
    } else {
        metadata$sample.id <- apply(metadata[, c(cell_type_var, sample_id_vars)], 1, paste, collapse = "_")
    }

    # Select the specified metadata columns and sample ID
    selected_metadata <- metadata[, c("sample.id", metadata_vars, cell_type_var)]

    if (split) {
        # Extract cell type information
        cell_types <- metadata[[cell_type_var]]

        # Split metadata by cell type
        split_metadata <- split(selected_metadata, cell_types)

        # Aggregate metadata by sample ID for each cell type
        aggregated_metadata_list <- lapply(split_metadata, function(df) {
            metadata.df <- df %>%
                group_by(sample.id) %>%
                summarise(across(everything(), ~ unique(.))) %>%
                column_to_rownames(var = "sample.id")

            # Replace any hyphens in rownames with "."
            rownames(metadata.df) <- gsub("-", ".", rownames(metadata.df))

            # return
            return(metadata.df)
        })

        # Return the list of aggregated metadata
        return(aggregated_metadata_list)
    } else {
        # Aggregate metadata by sample ID
        aggregated_metadata <- selected_metadata %>%
            group_by(sample.id) %>%
            summarise(across(everything(), ~ unique(.))) %>%
            column_to_rownames(var = "sample.id")

        # Replace any hyphens in rownames with "."
        rownames(aggregated_metadata) <- gsub("-", ".", rownames(aggregated_metadata))

        # Return the resulting data frame
        return(aggregated_metadata)
    }
}
