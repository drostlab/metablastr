#' @title Extract random loci from a set of genomes
#' @description In some cases, users may wish to extract sequences from randomly sampled loci of a particular length from a set of genomes.
#' This function allows users to specify a number of sequences
#' of a specified length that shall be randomly sampled from the genome.
#' The sampling rule is as follows:
#'For each locus independently sample:
#' \itemize{
#' \item 1) choose randomly (equal probability: see \code{\link{sample.int}} for details) from which of the given chromosomes the locus shall be sampled (\code{replace = TRUE}).
#' \item 2) choose randomly (equal probability: see \code{\link{sample.int}} for details) from which strand (plus or minus) the locus shall be sampled (\code{replace = TRUE}).
#' \item 3) randomly choose (equal probability: see \code{\link{sample.int}} the starting position of the locus in the sampled chromosome and strand (\code{replace = TRUE}).
#' }
#' @param sample_size a non-negative integer giving the number of loci that shall be sampled.
#' @param replace logical value indicating whether sampling should be with replacement. Default: \code{replace = TRUE}.
#' @param prob a vector of probability weights for obtaining the elements of the vector being sampled. Default is \code{prob = NULL}.
#' @param interval_width the length of the locus that shall be sampled.
#' @param subject_genomes a vector containing file paths to the reference genomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).
#' @param file_name name of the fasta file that stores the BLAST hit sequences. This name will only be used when \code{separated_by_genome = FALSE}.
#' @param separated_by_genome a logical value indicating whether or not hit sequences from different genomes should be stored in the same
#' output \code{fasta} file \code{separated_by_genome = FALSE} (default) or in separate \code{fasta} files \code{separated_by_genome = TRUE}.
#' @param update shall an existing \code{file_name} file be overwritten (\code{update = TRUE}; Default) or shall blast hit sequences be appended to the existing file (\code{update = FALSE})? 
#' @param path a folder path in which corresponding \code{fasta} output files shall be stored.
#' @author Hajk-Georg Drost
#' @export
extract_random_seqs_from_multiple_genomes <-
  function(sample_size,
           replace = TRUE, 
           prob = NULL,
           interval_width,
           subject_genomes,
           file_name = NULL,
           separated_by_genome = FALSE,
           update = TRUE,
           path = NULL) {
    
    if (any(!file.exists(subject_genomes)))
      stop(
        "At least one of the genome paths seems not to exist. Please check that all paths correspond to the correct location of the genome files.",
        call. = FALSE
      )
    
    if (!is.null(file_name) & separated_by_genome)
      stop(
        "A file name can only be specified when separated_by_genome = 'FALSE'. Otherwise separate species specific fasta files are being generated.",
        call. = FALSE
      )
    
    if (!is.null(path)){
      if (!file.exists(path)){
        message("The output folder '", path, "' does not seem to exist. A new folder will be created.")
        dir.create(path)
      }
    }
      
    message("Starting sequence extraction process ...")
    
    ifelse(
      separated_by_genome,
      seq_file_paths <-
        vector("character", 1),
      seq_file_paths <- vector("character", length(subject_genomes))
    )
    
    if (separated_by_genome) {
      for (i in seq_len(length(subject_genomes))) {
        # remove appendix *.fa from file name
        species_refined_name <-
          unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        message("Processing organism ", species_refined_name, " ...")
        
        extract_random_seqs_from_genome(
          size = sample_size,
          replace = replace,
          prob = prob,
          interval_width = interval_width,
          subject_genome = subject_genomes[i],
          append = FALSE,
          file_name = ifelse(is.null(path), paste0(species_refined_name[i],"_randomly_drawn_seqs_sample_size_", sample_size, "_interval_width_", interval_width, ".fa"), file.path(path, paste0(species_refined_name[i],"_randomly_drawn_seqs_sample_size_", sample_size, "_interval_width_", interval_width, ".fa")))
        )
        
      }
    } else {
      if (is.null(file_name)) {
        warning(
          "No file name was specified, thus the default file name '",paste0("randomly_drawn_seqs_sample_size_", sample_size, "_interval_width_", interval_width, ".fa"),"' was used.",
          call. = FALSE
        )
        file_name <- paste0("randomly_drawn_seqs_sample_size_", sample_size, "_interval_width_", interval_width, ".fa")
      }
      
      seq_file_paths <- ifelse(is.null(path),
                               file.path(getwd(), file_name),
                               file.path(path, file_name))
      
      if (update & file.exists(seq_file_paths))
        fs::file_delete(seq_file_paths)
      
      for (i in seq_len(length(subject_genomes))) {
         # remove appendix *.fa from file name
        species_refined_name <-
          unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        message("Processing organism ", species_refined_name, " ...")
        
        extract_random_seqs_from_genome(
          size = sample_size,
          replace = replace,
          prob = prob,
          interval_width = interval_width,
          subject_genome = subject_genomes[i],
          append = TRUE,
          file_name = seq_file_paths
        )
        
      }
    }
    message("Random sequence extraction process for all organisms finished without any problems. The result file was stored at '", seq_file_paths, "'.")
    return(seq_file_paths)
  }
