remove_multipart_transcripts <- function(x) {
  if (identical(x$strand, rev(x$strand))) {
    return(x)
  } else {
    return(x[1, ])
  }
}