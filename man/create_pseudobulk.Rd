% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_pseudobulk.R
\name{create_pseudobulk}
\alias{create_pseudobulk}
\title{Create Pseudo-bulked Assay data}
\usage{
create_pseudobulk(seurat_obj, cell_type_var, grouping_vars, split = FALSE)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{cell_type_var}{Chr name of metadata variable containing cell type annotation}

\item{grouping_vars}{Chr vector of metadata variables that define cell grouping strategy}

\item{split}{Lgl whether or not to split pseudo-bulk matrix by cell type}
}
\value{
Data frame, or list of data frames if \code{split == TRUE}
}
\description{
Create Pseudo-bulked Assay data
}
\examples{
# create_pseudobulk(
#   seurat_obj = Neurons,
#   cell_type_var = "CellType",
#   grouping_vars = c("Sex", "GEM"),
#   split = TRUE
# )
}
