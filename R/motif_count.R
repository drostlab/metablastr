#' @title Count the number of motifs in a set of sequences
#' @description Detect all matches of a given motif (typically a short string) in a (typically long sequence) reference sequence or set of reference sequences.
#' @param file a file path to the \code{fasta} file storing the sequences of interest.
#' @param motif a character string specifying the sequence motif that shall be searched.
#' @param max.mismatch the maximum number of mismatching letters allowed (see \code{\link[Biostrings]{matchPattern}} for details).
#' @param min.mismatch the minimum number of mismatching letters allowed (see \code{\link[Biostrings]{vcountPattern}} for details).
#' @param \dots additional arguments passed to \code{\link[Biostrings]{matchPattern}}.
#' @author Hajk-Georg Drost
#' @export
motif_count <- function(file, motif, max.mismatch = 0, min.mismatch = 0, ...) {
  
  if (!file.exists(file))
    stop("The file '",file,"' seems not to exist.", call. = FALSE)
  if (!is.character(motif))
    stop("The motif must be a character string.", call. = FALSE)
  if (length(motif) > 1)
    stop("Please specify only one motif at a time.", call. = FALSE)
    
  sequences <- Biostrings::readBStringSet(filepath = file, format = "fasta")
  count_motif <-
    Biostrings::vcountPattern(
      pattern = motif,
      subject = sequences,
      max.mismatch = max.mismatch,
      min.mismatch = min.mismatch,
      ...
    )
  names(count_motif) <- names(sequences)
  return(count_motif)
}

