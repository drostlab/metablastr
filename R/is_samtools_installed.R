is_samtools_installed <- function() {
  # test if a valid samtools version is installed
  tryCatch({
    sys_out <-
      system("samtools help", intern = TRUE)
  }, error = function(e)
    stop(
      "It seems like you don't have samtools installed locally on your machine or the PATH variable to the samtools program is not set correctly. Please install samtools (http://www.htslib.org/download/) before running this command.",
      call. = FALSE
    ))
  
  if (any(stringr::str_detect(sys_out[2], "samtools")))
    return(TRUE)
  
}
