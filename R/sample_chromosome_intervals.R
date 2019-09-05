#' @title Helper function to sample random intervals of length \code{interval_width} from chromosomes
#' @description This function samples a random locus given the chromosome lenght (\code{chr_size}) 
#' and \code{interval_width} from a chromosome
#' @param chr_size an integer value defining the length of the chromosome.
#' @param interval_width the length of the interval that shall be sampled.
#' @param strand from which strand should the interval be sampled? Options
#' are:
#' \itemize{
#' \item \code{strand = "plus"}
#' \item \code{strand = "minus"}
#' }
#' @param size a non-negative integer giving the number of sequences that shall be sampled from random genomic loci.
#' @param replace logical value indicating whether sampling should be with replacement. Default: \code{replace = TRUE}.
#' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Default is \code{prob = NULL}.
#' @author Hajk-Georg Drost
sample_chromosome_intervals <-
  function(chr_size,
           interval_width,
           strand,
           size,
           replace = TRUE,
           prob = NULL) {
    
  if (!is.element(strand, c("plus", "minus")))
    stop("The 'strand' argument can only be specified as strand = 'plus' or strand = 'minus'.", call. = FALSE)
  
  if (chr_size - interval_width < 0) {
    message("It seems that the chromosome length is smaller than the sampled locus. Thus, this chromosome will be omitted.")
    res <-
      tibble::tibble(start = NA,
                     end = NA,
                     width = NA,
                     strand = NA)
    return(res)
    }
  
  # random start position in chromosome
    random_start <-
      sample.int(ifelse((chr_size - interval_width) == 0, 1, chr_size - interval_width),
                 size,
                 replace = replace,
                 prob = prob)
  # compute random end position in chromosome given interval width
  random_end <- random_start + interval_width - 1
  
  start <- end <- NULL
  if (strand == "plus") {
    res <-
      tibble::tibble(start = random_start,
                     end = random_end,
                     width = end - start + 1,
                     strand = rep("plus", length(random_start)))
  }
  
  if (strand == "minus") {
    res <-
      tibble::tibble(start = random_end,
                     end = random_start,
                     width = start - end + 1,
                     strand = rep("minus", length(random_start)))
  }  
  return(res)
}




