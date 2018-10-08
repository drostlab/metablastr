#' @title Retrieve the best hits across species from a BLAST table
#' @description This function aims to retrieve the best blast hits for each query sequence
#' by filtering a \code{blast_tbl} using the following criteria.
#' A best hit is defined as (fulfilling all three critaria):
#' \itemize{
#' \item max(alig_length): only the hit having the longest alignment length is retained.
#' \item \code{qcovhsp >= min_qcovhsp}: only hits that have a query coverage of at least \code{min_qcovhsp} are retained.
#' \item max(bit_score): only the hit having the highest bit-score is retained.  
#' }  
#' @param blast_tbl a BLAST table generated with \code{\link{detect_homologs_proteome_to_proteome}} or \code{\link{detect_homologs_cds_to_cds}}.
#' @param min_qcovhsp minimum query coverage of the hit in percent \code{10..100} that shall be retained. Default value is set to \code{min_qcovhsp = 50} (= a best hit alignment must have at least 50% query coverage).
#' @author Hajk-Georg Drost
#' @export

filter_best_hits <- function(blast_tbl, min_qcovhsp = 50) {
  
  if (!dplyr::between(min_qcovhsp, 10, 100))
    stop("Please provide a min_qcovhsp value between 10 and 100.", call. = FALSE)
  
  if (nrow(blast_tbl) == 0)
    stop("Please provide a blast_tbl that contains at least one row.", call. = FALSE)
  
  alig_length <- qcovhsp <- bit_score <- species <- query_id <- NULL
  
  message("Retrieving best blast hits using the following criteria: ")
  message(" 1) the query coverage ('qcovhsp') of the hit must be at least greater than ", qcovhsp)
  message(" 2) select the blast hit with maximum 'alig_length'")
  message(" 3) select the blast hit that in addition has the maximum bit_score")
  message("--------")
  message("Number of hits before filtering: ", nrow(blast_tbl))
  blast_tbl <- dplyr::filter(blast_tbl, qcovhsp >= min_qcovhsp)

  filter_best_hits <- function(x) {
    res <-
      dplyr::filter(x, max(alig_length), max(bit_score))
    return(res)
  }
  
  best_hit_df <- dplyr::do(dplyr::group_by(blast_tbl, species, query_id), filter_best_hits(.))
  
  message("Number of best hits after filtering: ", nrow(best_hit_df))
  
  if (nrow(best_hit_df) == 0)
    stop("The filter process resultet in 0 best hits. Please provide more liberal filter criteria to retrieve a best hit table.", call. = FALSE)
  
  return(best_hit_df)
}

