#' @title Test motif enrichment in a set of non-random versus random sequences
#' @description Compare the number of motifs in a set of non-random versus random sequences.
#' The resulting values are then tested for enrichment of certain motifs in real sequences compared to random sequences. 
#' Several tests statistics and approaches are available to quantify significant motif enrichment.
#' @param real_seqs a file path to the \code{fasta} file storing the non-random set of sequences.
#' @param random_seqs a file path to the \code{fasta} file storing the random set of sequences, e.g. generated with \code{\link{extract_random_seqs_from_genome}}.
#' @param motifs a character vector storing a set of motifs that shall be counted within respective sequences.
#' @param test tests statistics and models to quantify significant motif enrichment. Options are:
#' \itemize{
#' \item \code{test = "fisher"}: Fisher's Exact Test for Count Data (see \code{link[stats]{fisher.test}} for details).
#' }
#' @param alternative indicates the alternative hypothesis and must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just the initial letter. Only used in the 2 by 2 case.
#' @param max.mismatch the maximum number of mismatching letters allowed (see \code{\link[Biostrings]{matchPattern}} for details).
#' @param min.mismatch the minimum number of mismatching letters allowed (see \code{\link[Biostrings]{vcountPattern}} for details).
#' @param \dots additional arguments passed to \code{\link[Biostrings]{matchPattern}}.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{motif_count}}, \code{\link{extract_random_seqs_from_genome}}, \code{\link{extract_hit_seqs_from_genomes}}
#' @export

motif_enrichment <-
  function(real_seqs,
           random_seqs,
           motifs,
           test = "fisher",
           alternative = "less",
           max.mismatch = 0,
           min.mismatch = 0,
           ...) {
    if (!is.element(test, c("fisher")))
      stop(
        "The test type '",
        test,
        "' is not available for this function. Please consult the documentation for details.",
        call. = FALSE
      )
    
    '.' <- NULL
    
    motif_compare_tbl <-
      motif_compare(real_seqs = real_seqs,
                    random_seqs = random_seqs,
                    motifs = motifs,
                    max.mismatch = max.mismatch,
                    min.mismatch = min.mismatch,
                    ...
                    )
    
    if (test == "fisher") {
      fisher_p_val <- function(x) {
        seqtype <- NULL
        real_seqs_motif_count <-
          dplyr::filter(x, seqtype == "real")$motif_count
        real_seqs_seq_count <-
          dplyr::filter(x, seqtype == "real")$seq_count
        random_seqs_motif_count <-
          dplyr::filter(x, seqtype == "random")$motif_count
        random_seqs_seq_count <-
          dplyr::filter(x, seqtype == "random")$seq_count
        
        count_tbl <-
          matrix(
            c(
              real_seqs_seq_count,
              real_seqs_motif_count,
              random_seqs_seq_count,
              random_seqs_motif_count
            ),
            nrow = 2,
            dimnames = list(
              c("real_seq_count", "real_motif_count"),
              c("random_seq_count", "random_motif_count")
            )
          )
        
        return(stats::fisher.test(count_tbl, alternative = alternative)$p.value)
      }
      
      motif <- NULL
      res <-
        dplyr::do(
          dplyr::group_by(motif_compare_tbl, motif),
          data.frame(fisher_pval = fisher_p_val(.))
        )
    }
    return(res)
  }



