#' @title Retrieve a core set of best hits that is shared across at least N species
#' @description This function aims to retrieve a core set of best blast hits for each query sequence that is shared across all species. In other words, only query sequences
#' that generated blast hits in all species in the input \code{blast_tbl} were retained.
#' by filtering a \code{blast_tbl} using the following criteria.
#' A best hit is defined as (fulfilling all three critaria):
#' \itemize{
#' \item max(alig_length): only the hit having the longest alignment length is retained.
#' \item \code{qcovhsp >= min_qcovhsp}: only hits that have a query coverage of at least \code{min_qcovhsp} are retained.
#' \item max(bit_score): only the hit having the highest bit-score is retained.
#' }
#' @param blast_tbl a BLAST table generated with \code{\link{detect_homologs_proteome_to_proteome}} or \code{\link{detect_homologs_cds_to_cds}}.
#' @param min_qcovhsp minimum query coverage of the hit in percent \code{0..100} that shall be retained. Default value is set to \code{min_qcovhsp = 50} (= a best hit alignment must have at least 50 percent query coverage).
#' @param n_species minimum number of species in which a query sequence needs to generates valid blast hits. 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{filter_best_hits}}, \code{\link{filter_homologs_core_set}}
#' @export

filter_homologs_core_set_n_species <- function(blast_tbl, min_qcovhsp = 50, n_species) {
  
  message("Retrieve core set of homologs (= query hits shared across ",n_species," species).")
  
  qcovhsp <- query_id <- '.' <- NULL
  
  blast_tbl_processed <- dplyr::filter(blast_tbl, qcovhsp >= min_qcovhsp)
  
  # total number of species in the BLAST table
  total_n_species <-
    length(names(table(unique(
      blast_tbl$species
    ))))
  
  if (length(total_n_species) == 0 || is.na(total_n_species) || n_species > total_n_species)
    stop(
      "It seems like your blast_tbl does not contain any species information or you specified n_species larger than uniqe species available in blast_tbl.",
      call. = FALSE
    )
  
  # create empty row
  col_names_df <- names(blast_tbl)
  tmp_NULL_df <- vector("list", length(col_names_df))
  names(tmp_NULL_df) <- col_names_df
  tmp_NULL_df <- lapply(tmp_NULL_df, function(x)
    x <- NA)
  empty_row_df <- dplyr::bind_cols(tmp_NULL_df)
  
  
  filter_species_hits <- function(x) {
    nspecies <- length(names(table(unique(as.character(x$species)))))
    if (nspecies >= n_species) {
      return(x)
    } else {
      return(empty_row_df)
    }
  }
  
  query_id <- gene_name <- NULL
  blast_tbl_processed_hit_in_all_species <-
    dplyr::do(dplyr::group_by(blast_tbl_processed, gene_name),
              filter_species_hits(.))
  
  blast_tbl_processed_hit_in_all_species <-
    dplyr::filter(blast_tbl_processed_hit_in_all_species, !is.na(gene_name))
  
  message("Number of core best hits shared across all species after filtering: ", nrow(blast_tbl_processed_hit_in_all_species))
  
  if (nrow(blast_tbl_processed_hit_in_all_species) == 0) {
    stop(
      "The filter process resultet in 0 best hits. Please provide more liberal filter criteria to retrieve a core set of hits shared across all species.",
      call. = FALSE
    )
  }
  
  return(blast_tbl_processed_hit_in_all_species)
}
