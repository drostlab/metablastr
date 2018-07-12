#' @title 
#' @description 
#' @param size a non-negative integer giving the number of loci that shall be sampled.
#' @param interval_width
#' @param subject_genome
#' @param file_name
#' @author 
#' @export

extract_random_seqs_from_genome <-
  function(size,
           interval_width,
           subject_genome,
           file_name = NULL) {
    
    if (any(!file.exists(subject_genome)))
      stop(
        "The genome path seems not to exist. Please check that the path corresponds to the correct location of the genome file.",
        call. = FALSE
      )
    
      if(is.null(file_name)) {
        warning("No file name was specified, thus the default file name '", paste0(basename(subject_genome), "_random_seqs.fa"), "' was used.", call. = FALSE)
        file_name <- file.path(getwd(), paste0(basename(subject_genome), "_random_seqs.fa"))
      }
      
    strand <- chr <- NULL
    
        # remove appendix *.fa from file name
        species_refined_name <- unlist(stringr::str_split(basename(subject_genome), "[.]"))[1]
        message("Processing organism ", species_refined_name, " ...")
        imported_genome_i <- biomartr::read_genome(subject_genome)
        
        # only retain chromosome names that are present in both: genome and BLAST table
        chr_names <-
          unlist(lapply(names(imported_genome_i), function(x)
            unlist(stringr::str_split(x, " ")[1])))
        if (length(chr_names) == 0)
          stop("No chromosomes were found in this genome.", call. = FALSE)
        
        names(imported_genome_i) <- chr_names
        
        res <- vector("list", length(size))
        
        for (i in seq_len(size)) {
          sample_chromsome <- sample.int(length(chr_names), 1, replace = TRUE)
          sample_strand <- sample.int(2, 1, replace = TRUE, prob = c(0.5, 0.5))
          
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
            sample_i <- dplyr::mutate(sample_i, chr = chr_names[sample_chromsome])
            res[i] <- list(sample_i)
          }
        }
        
        random_coordinates <- dplyr::bind_rows(res)
        
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