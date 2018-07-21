#' @title Count the number of motifs in a set of sequences and generate a motif count matrix
#' @description Detect all matches of a given motif (typically a short string) in a (typically long sequence) reference sequence or set of reference sequences. For each sequence and each motif a count matrix is generated and returned.
#' @param file a file path to the \code{fasta} file storing the sequences of interest.
#' @param motifs a character string specifying the sequence motifs that shall be screened and counted.
#' @param max.mismatch the maximum number of mismatching letters allowed (see \code{\link[Biostrings]{matchPattern}} for details).
#' @param min.mismatch the minimum number of mismatching letters allowed (see \code{\link[Biostrings]{vcountPattern}} for details).
#' @param \dots additional arguments passed to \code{\link[Biostrings]{matchPattern}}.
#' @author Hajk-Georg Drost
#' @export
motif_count_matrix <- function(file, motifs, max.mismatch = 0, min.mismatch = 0, ...) {
  
  if (!file.exists(file))
    stop("Please provide a valid path to your input fasta file.", call. = FALSE)
  
  input_file <- Biostrings::readBStringSet(file)
  
  res <- matrix(NA_integer_, length(motifs), length(input_file@ranges@NAMES))
  message("Generating motif count matrix ...")
  for (i in seq_len(length(motifs))) {
    message(" -> processing motif '", motifs[i], "'")
    
    motif_count_i <- metablastr::motif_count(
      file = file,
      motif = motifs[i],
      max.mismatch = max.mismatch,
      min.mismatch = min.mismatch,
      ...
    )
    res[i, ] <- motif_count_i
  }
  res <- as.data.frame(res)
  names(res) <- input_file@ranges@NAMES
  res <- cbind(motifs, res)
  names(res)[1] <- "motifs"
  message("Count matrix successfully generated!")
  return(tibble::as.tibble(res))
}