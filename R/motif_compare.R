#' @title Count the number of motifs in a set of non-random versus random sequences 
#' @description Compare the number of motifs in a set of non-random versus random sequences.
#' The resulting values can then be used to test the enrichment of certain motifs in real sequences compared to random sequences.
#' @param real_seqs a file path to the \code{fasta} file storing the non-random set of sequences.
#' @param random_seqs a file path to the \code{fasta} file storing the random set of sequences, e.g. generated with \code{\link{extract_random_seqs_from_genome}}.
#' @param motifs a character vector storing a set of motifs that shall be counted within respective sequences.
#' @param max.mismatch the maximum number of mismatching letters allowed (see \code{\link[Biostrings]{matchPattern}} for details).
#' @param min.mismatch the minimum number of mismatching letters allowed (see \code{\link[Biostrings]{vcountPattern}} for details).
#' @param \dots additional arguments passed to \code{\link[Biostrings]{matchPattern}}.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{motif_count}}, \code{\link{extract_random_seqs_from_genome}}, \code{\link{extract_hit_seqs_from_genomes}}
#' @export

motif_compare <- function(real_seqs, random_seqs, motifs, max.mismatch = 0, min.mismatch = 0, ...) {
  
  res <- vector("list", length(motifs))
  
  for(i in seq_len(length(motifs))) {
    
    real_seqs_motifcount <- sum(motif_count(file = real_seqs, motif = motifs[i]), max.mismatch = max.mismatch, min.mismatch = min.mismatch, ...)
    real_seqs_import <- Biostrings::readBStringSet(real_seqs)
    real_seqs_seqcount <- length(real_seqs_import)
    random_seqs_motifcount <- sum(motif_count(file = random_seqs, motif = motifs[i]), max.mismatch = max.mismatch, min.mismatch = min.mismatch, ...)
    random_seqs_import <- Biostrings::readBStringSet(random_seqs)
    random_seqs_seqcount <- length(random_seqs_import)
    
    res[i] <- list(tibble::tibble(
      motif = c(motifs[i], motifs[i]),
      seqtype = c("real", "random"),
      motif_count = c(real_seqs_motifcount, random_seqs_motifcount),
      seq_count = c(real_seqs_seqcount, random_seqs_seqcount)
    ))
  }
  
  return(dplyr::bind_rows(res))
}