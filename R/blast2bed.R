#' @title Export BLAST output table to BED file
#' @description Taking a \code{blast_tbl} as input this function exports
#' the following columns as BED file.
#' \itemize{
#' \item \code{chromosome} is defined by the \code{blast_tbl} column \code{subject_id}.
#' \item \code{start} is defined by the \code{blast_tbl} column \code{s_start}.
#' \item \code{end} is defined by the \code{blast_tbl} column \code{s_end}.
#' }
#' @param blast_tbl a BLAST table generated with any \code{blast_*()} function.
#' @param output file path or name of the output file.
#' @author Hajk-Georg Drost
#' @export
blast2bed <- function(blast_tbl, output) {
  
  res <- dplyr::data_frame(chromosome = blast_tbl$subject_id, 
                           start      = blast_tbl$s_start, 
                           end        = blast_tbl$s_end)
  
  readr::write_tsv(res, output, col_names = FALSE)
}