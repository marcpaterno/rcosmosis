#' @importFrom methods setOldClass
#' @exportClass cosmosis_tbl
setOldClass(c("cosmosis_tbl", "tbl_df", "tbl", "data.frame"))

#' `cosmosis_tbl` class
#'
#' @description
#' The `cosmosis_tbl` class is a subclass of [`tbl_df`][tibble::tbl_df()],
#' created in order to have different default behaviour.
#'
#' `cosmosis_tbl` provides functions for a few generic methods:
#'     summary: provides point estimates and credible regions for each
#'              sampled parameter
#'
#' @section Properties of `cosmosis_tbl`:
#'
#' Objects of class `tbl_df` have:
#' * A `class` attribute of `c("cosmosis_tbl", tbl_df", "tbl", "data.frame")`.
#'
#' @name cosmosis_tbl-class
NULL
