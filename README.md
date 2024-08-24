
<!-- README.md is generated from README.Rmd. Please edit that file -->

# miSCutils

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/gtwa-bio/miSCutils)](https://github.com/gtwa-bio/miSCutils/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/gtwa-bio/miSCutils)](https://github.com/gtwa-bio/miSCutils/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/gtwa-bio/miSCutils/graph/badge.svg)](https://app.codecov.io/gh/gtwa-bio/miSCutils)
<!-- badges: end -->

The goal of `miSCutils` is to collect the miscellaneous functions that
tend to get written and reused during single-cell analysis. Rather than
copy/pasting and using slightly modified versions repeatedly, this
package will serve as a common collection for version controlling and
re-distributing these utilities elsewhere. It **will not** contain any
novel data analysis or statistical methods. Functions will primarily be
wrappers for established packages and workflows.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `miSCutils` from
[GitHub](https://github.com/gtwa-bio/PMIDlotteRy) with:

``` r
install.packages("devtools")
devtools::install_github("gtwa-bio/PMIDlotteRy")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## !! WIP !! ##
```

## Code of Conduct

Please note that the `miSCutils` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.19/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://gtwa-bio.github.io/miSCutils) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.19/biocthis)*.

## Acknowledgements

Please note that the `miSCutils` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.
