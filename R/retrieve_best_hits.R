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
#' @param min_qcovhsp minimum query coverage of the hit in percent \code{0..100} that shall be retained. Default value is set to \code{min_qcovhsp = 50} (= a best hit alignment must have at least 50% query coverage).
#' @author Hajk-Georg Drost
#' @export

retrieve_best_hits <- function(blast_tbl, min_qcovhsp = 50) {
  
  alig_length <- qcovhsp <- bit_score <- species <- query_id <- NULL
  filter_best_hits <- function(x) {
    res <-
      dplyr::filter(x, max(alig_length), qcovhsp >= min_qcovhsp, max(bit_score))
    return(res)
  }
  
  best_hit_df <- dplyr::do(dplyr::group_by(blast_tbl, species, query_id), filter_best_hits(.))
  return(best_hit_df)
}

