#' @title Extract random loci from a genome of interest 
#' @description This function allows users to specify a number of sequences
#' of a specified length that shall be randomly sampled from the genome.
#' The sampling rule is as follows:
#'For each locus independently sample:
#' \itemize{
#' \item 1) choose randomly (equal probability: see \code{\link{sample.int}} for details) from which of the given chromosomes the locus shall be sampled (\code{replace = TRUE}).
#' \item 2) choose randomly (equal probability: see \code{\link{sample.int}} for details) from which strand (plus or minus) the locus shall be sampled (\code{replace = TRUE}).
#' \item 3) randomly choose (equal probability: see \code{\link{sample.int}} the starting position of the locus in the sampled chromosome and strand (\code{replace = TRUE}).
#' }
#' @param size a non-negative integer giving the number of loci that shall be sampled.
#' @param replace logical value indicating whether sampling should be with replacement. Default: \code{replace = TRUE}.
#' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Default is \code{prob = NULL}.
#' @param interval_width the length of the locus that shall be sampled.
#' @param subject_genome file path to the \code{fasta} file storing the subject genome.
#' @param file_name a name of the output \code{fasta} file that will store the sequences of the randomly
#' sampled loci.
#' @param append shall new random sequences be added to an existing \code{file_name} (\code{append = TRUE})
#'  or should an exosting \code{file_name} be removed before storing new random sequences (\code{append = FALSE}; Default)?  
#' @author Hajk-Georg Drost
#' @export

extract_random_seqs_from_genome <-
  function(size,
           replace = TRUE,
           prob = NULL,
           interval_width,
           subject_genome,
           file_name = NULL,
           append = FALSE) {
    
    if (!file.exists(subject_genome))
      stop(
        "The genome path seems not to exist. Please check that the path corresponds to the correct location of the genome file.",
        call. = FALSE
      )
    
    if (length(subject_genome) > 1)
      stop("Please provide only one subject genome.", call. = FALSE)
    
      if(is.null(file_name)) {
        warning("No file name was specified, thus the default file name '", paste0(basename(subject_genome), "_random_seqs.fa"), "' was used.", call. = FALSE)
        file_name <- file.path(getwd(), paste0(basename(subject_genome), "_random_seqs.fa"))
      }
      
    if (length(interval_width) > 1)
      stop("Please provide only one width.", call. = FALSE)
    
    if (file.exists(file_name) & !append)
      file.remove(file_name)
    
    strand <- chr <- NULL
    
    # remove appendix *.fa from file name
    split_name <- unlist(stringr::str_split(basename(subject_genome), "[.]"))
    if (length(split_name) > 1) {
      species_refined_name <- paste0(split_name[-length(split_name)], collapse = ".")
    } else {
      species_refined_name <- split_name[1]
    }
    
    message("Processing organism ", species_refined_name, " ...")
        imported_genome_i <- biomartr::read_genome(subject_genome)
        
        # only retain chromosome names that are present in both: genome and BLAST table
        chr_names <-
          unlist(lapply(names(imported_genome_i), function(x) {
            new_name <- unlist(stringr::str_split(x, " "))[1]
            new_name <- unlist(stringr::str_replace_all(new_name, "[.]", "_"))
            new_name <- stringr::str_trim(new_name, side = "both")
            return(new_name)
          }))
            
        if (length(chr_names) == 0)
          stop("No chromosomes were found in this genome.", call. = FALSE)
        
        names(imported_genome_i) <- chr_names
        
        res <- vector("list", length(size))
        
        for (i in seq_len(size)) {
          sample_chromsome <- sample.int(length(chr_names), 1, replace = replace, prob = prob)
          sample_strand <- sample.int(2, 1, replace = replace, prob = c(0.5, 0.5))
          
          if (sample_strand == 1){
            sample_i <- sample_chromosome_intervals(
              chr_size = imported_genome_i[chr_names[sample_chromsome]]@ranges@width,
              interval_width = interval_width,
              strand = "plus",
              size = 1
            )
            
            sample_i <- dplyr::mutate(sample_i, chr = chr_names[sample_chromsome])
            res[i] <- list(sample_i)
          }
          
          if (sample_strand == 2){
            sample_i <- sample_chromosome_intervals(
              chr_size = imported_genome_i[chr_names[sample_chromsome]]@ranges@width,
              interval_width = interval_width,
              strand = "minus",
              size = 1
            )
            
                sample_i <-
                  dplyr::mutate(sample_i, chr = chr_names[sample_chromsome])
                res[i] <- list(sample_i)
          }
        }
        
        random_coordinates <- dplyr::bind_rows(res)
        random_coordinates <- stats::na.omit(random_coordinates)
        
          # for each chromosome separately
          for (j in seq_len(length(chr_names))) {
            
            # divide by plus and minus strand
            species_specific_random_plus_strand <-
              dplyr::filter(random_coordinates, chr == chr_names[j], strand == "plus")
            species_specific_random_minus_strand <-
              dplyr::filter(random_coordinates, chr == chr_names[j], strand == "minus")
            
            if (nrow(species_specific_random_plus_strand) > 0) {
              message("  -> Extracting random sequences from the plus strand of ", chr_names[j], " ...")
              
              seqs_plus <- Biostrings::extractAt(
                imported_genome_i[chr_names[j]],
                at = IRanges::IRanges(
                  start = species_specific_random_plus_strand$start,
                  end = species_specific_random_plus_strand$end,
                  names = paste0(
                    species_refined_name,
                    "_",
                    chr_names[j],
                    "_",
                    species_specific_random_plus_strand$start,
                    "_",
                    species_specific_random_plus_strand$end,
                    "_strand_plus"
                  )
                )
              )
              
              Biostrings::writeXStringSet(
                seqs_plus@unlistData,
                filepath = file_name,
                format = "fasta",
                append = TRUE
              )
              
            } else {
              message("  -> No random sequences were chosen from the plus strand of ", chr_names[j])
            }
            
            if (nrow(species_specific_random_minus_strand) > 0) {
              message("  -> Extracting random sequences from the minus strand of ", chr_names[j], " ...")
              
              seqs_minus <- Biostrings::extractAt(
                Biostrings::reverse(imported_genome_i[chr_names[j]]),
                at = IRanges::IRanges(
                  start = imported_genome_i[chr_names[j]]@ranges@width - species_specific_random_minus_strand$start + 1,
                  end = imported_genome_i[chr_names[j]]@ranges@width - species_specific_random_minus_strand$end + 1,
                  names = paste0(
                    species_refined_name,
                    "_",
                    chr_names[j],
                    "_",
                    species_specific_random_minus_strand$start,
                    "_",
                    species_specific_random_minus_strand$end,
                    "_strand_minus"
                  )
                )
              )
              
              Biostrings::writeXStringSet(
                seqs_minus@unlistData,
                filepath = file_name,
                format = "fasta",
                append = TRUE
              )
              
            } else {
              message("  -> No random sequences were chosen from the minus strand of ", chr_names[j])
            }
          }
        
    message("Sequence extraction process of random loci finished without any problems.")
    return(file_name)
  }
