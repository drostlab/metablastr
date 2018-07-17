#' @title Extract sequences of BLAST hits in respective genomes and store it as \code{fasta} file(s)
#' @description In some cases, users may wish to extract the sequences of the respective blast hit
#' stored in a \code{blast_tbl}. This function enables to quickly extract such sequences and store
#' them in one general or multiple genome specific \code{fasta} file(s).
#' @param blast_tbl a BLAST table generated with \code{\link{blast_genomes}}.
#' @param subject_genomes a vector containing file paths to the reference genomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).
#' @param file_name name of the fasta file that stores the BLAST hit sequences. This name will only be used when \code{separated_by_genome = FALSE}.
#' @param separated_by_genome a logical value indicating whether or not hit sequences from different genomes should be stored in the same
#' output \code{fasta} file \code{separated_by_genome = FALSE} (default) or in separate \code{fasta} files \code{separated_by_genome = TRUE}.
#' @param path a folder path in which corresponding \code{fasta} output files shall be stored.
#' @author Hajk-Georg Drost
#' @export
extract_hit_seqs_from_genomes <-
  function(blast_tbl,
           subject_genomes,
           file_name = NULL,
           separated_by_genome = FALSE,
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
    
    message("Starting sequence extraction process ...")
    
    ifelse(
      separated_by_genome,
      seq_file_paths <-
        vector("character", 1),
      seq_file_paths <- vector("character", length(subject_genomes))
    )
    
    species <- s_strand <- subject_id <- NULL
    
    available_species <- names(table(blast_tbl$species))
    
    if (separated_by_genome) {
      for (i in seq_len(length(subject_genomes))) {
        # remove appendix *.fa from file name
        species_refined_name <-
          unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        message("Processing organism ", species_refined_name, " ...")
        if (!is.element(species_refined_name, available_species)) {
          message("Organism ",
                  species_refined_name,
                  " didn't have any BLAST hits.")
        } else {
          imported_genome_i <- biomartr::read_genome(subject_genomes[i])
          
          species_specific_blast_tbl <-
            dplyr::filter(blast_tbl, species == species_refined_name)
          
          if (nrow(species_specific_blast_tbl) > 0) {
            strand <-
              ifelse(
                species_specific_blast_tbl$s_end - species_specific_blast_tbl$s_start + 1 > 0,
                "plus",
                "minus"
              )
            
            species_specific_blast_tbl <-
              dplyr::mutate(species_specific_blast_tbl, s_strand = strand)
            
            imported_genome_i_names <-
              unlist(lapply(stringr::str_trim(names(imported_genome_i), side = "both"), function(x)
                unlist(stringr::str_split(x, " ")[1])))
            
            print(imported_genome_i_names)
            print(names(table(
              species_specific_blast_tbl$subject_id
            )))
            
            
            # only retain chromosome names that are present in both: genome and BLAST table
            chr_names <-
              dplyr::intersect(imported_genome_i_names, names(table(
                species_specific_blast_tbl$subject_id
              )))
            
            if (length(chr_names) == 0)
              stop(
                "It seems that the chromosome names used in the input 'blast_tbl' and in the input genome fasta file do not match. Please make sure that the chromosome names match in both cases.",
                call. = FALSE
              )
            
            # construct output file paths
            seq_file_paths[i] <- ifelse(is.null(path),
                                        file.path(
                                          getwd(),
                                          paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                                        ),
                                        file.path(
                                          path,
                                          paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                                        ))
            
            # for each chromosome separately
            for (j in seq_len(length(chr_names))) {
              message("  -> Extracting BLAST hit sequences from ",
                      chr_names[j],
                      " ...")
              # divide by plus and minus strand
              
              species_specific_blast_tbl_plus_strand <-
                dplyr::filter(species_specific_blast_tbl,
                              subject_id == chr_names[j],
                              s_strand == "plus")
              species_specific_blast_tbl_minus_strand <-
                dplyr::filter(
                  species_specific_blast_tbl,
                  subject_id == chr_names[j],
                  s_strand == "minus"
                )
              
              if (!identical(
                unique(species_specific_blast_tbl_plus_strand$s_len),
                imported_genome_i[chr_names[j]]@ranges@width
              ))
                stop(
                  "The chromosome length of ",
                  chr_names[j],
                  " (",
                  imported_genome_i[chr_names[j]]@ranges@width,
                  ") does not match the chromosome length s_len in the blast_tbl (",
                  unique(species_specific_blast_tbl_plus_strand$s_len),
                  "). Thus, BLAST hit sequences cannot be extracted due to indexing shifts...",
                  call. = FALSE
                  
                )
              
              seqs_plus <- Biostrings::extractAt(
                imported_genome_i[chr_names[j]],
                at = IRanges::IRanges(
                  start = species_specific_blast_tbl_plus_strand$s_start,
                  end = species_specific_blast_tbl_plus_strand$s_end,
                  names = paste0(
                    species_specific_blast_tbl_plus_strand$species,
                    "_",
                    species_specific_blast_tbl_plus_strand$subject_id,
                    "_",
                    species_specific_blast_tbl_plus_strand$s_start,
                    "_",
                    species_specific_blast_tbl_plus_strand$s_end,
                    "_strand_plus"
                  )
                )
              )
              
              seqs_minus <- Biostrings::extractAt(
                Biostrings::reverse(imported_genome_i[chr_names[j]]),
                at = IRanges::IRanges(
                  start = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_start + 1,
                  end = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_end + 1,
                  names = paste0(
                    species_specific_blast_tbl_minus_strand$species,
                    "_",
                    species_specific_blast_tbl_minus_strand$subject_id,
                    "_",
                    species_specific_blast_tbl_minus_strand$s_start,
                    "_",
                    species_specific_blast_tbl_minus_strand$s_end,
                    "_strand_minus"
                  )
                )
              )
              
              # store sequences in fasta files
              Biostrings::writeXStringSet(
                seqs_plus@unlistData,
                filepath = ifelse(
                  is.null(path),
                  file.path(
                    getwd(),
                    paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                  ),
                  file.path(
                    path,
                    paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                  )
                ),
                format = "fasta",
                append = TRUE
              )
              
              Biostrings::writeXStringSet(
                seqs_minus@unlistData,
                filepath = ifelse(
                  is.null(path),
                  file.path(
                    getwd(),
                    paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                  ),
                  file.path(
                    path,
                    paste0(subject_genomes[i], "_blast_tbl_sequences.fa")
                  )
                ),
                format = "fasta",
                append = TRUE
              )
            }
          }
        }
      } 
    } else {
      if (is.null(file_name)) {
        warning(
          "No file name was specified, thus the default file name 'blast_hit_seqs.fa' was used.",
          call. = FALSE
        )
        file_name <- "blast_hit_seqs.fa"
      }
      
      seq_file_paths <- ifelse(is.null(path),
                               file.path(getwd(), file_name),
                               file.path(path, file_name))
      
      for (i in seq_len(length(subject_genomes))) {
        message("Processing organism ",
                basename(subject_genomes[i]),
                " ...")
        # remove appendix *.fa from file name
        species_refined_name <-
          unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        if (!is.element(species_refined_name, available_species)) {
          message("Organism ",
                  species_refined_name,
                  " didn't have any BLAST hits.")
        } else {
          imported_genome_i <- biomartr::read_genome(subject_genomes[i])
          
          
          species_specific_blast_tbl <-
            dplyr::filter(blast_tbl, species == species_refined_name)
          
          if (nrow(species_specific_blast_tbl) > 0) {
            strand <-
              ifelse(
                species_specific_blast_tbl$s_end - species_specific_blast_tbl$s_start + 1 > 0,
                "plus",
                "minus"
              )
            species_specific_blast_tbl <-
              dplyr::mutate(species_specific_blast_tbl, s_strand = strand)
            
        
            imported_genome_i_names <-
              unlist(lapply(stringr::str_trim(names(imported_genome_i), side = "both"), function(x)
                unlist(stringr::str_split(x, " ")[1])))
            
            print(imported_genome_i_names)
            print(names(table(
              species_specific_blast_tbl$subject_id
            )))
            
            # only retain chromosome names that are present in both: genome and BLAST table
            chr_names <-
              dplyr::intersect(names(imported_genome_i), names(table(
                species_specific_blast_tbl$subject_id
              )))
            
            if (length(chr_names) == 0)
              stop(
                "It seems that the chromosome names used in the input 'blast_tbl' and in the input genome fasta file do not match. Please make sure that the chromosome names match in both cases.",
                call. = FALSE
              )
            
            # for each chromosome separately
            for (j in seq_len(length(chr_names))) {
              message("  -> Extracting BLAST hit sequences from ",
                      chr_names[j],
                      " ...")
              
              # divide by plus and minus strand
              species_specific_blast_tbl_plus_strand <-
                dplyr::filter(species_specific_blast_tbl,
                              subject_id == chr_names[j],
                              s_strand == "plus")
              species_specific_blast_tbl_minus_strand <-
                dplyr::filter(
                  species_specific_blast_tbl,
                  subject_id == chr_names[j],
                  s_strand == "minus"
                )
              
              if (!identical(
                unique(species_specific_blast_tbl_plus_strand$s_len),
                imported_genome_i[chr_names[j]]@ranges@width
              ))
                stop(
                  "The chromosome length of ",
                  chr_names[j],
                  " (",
                  imported_genome_i[chr_names[j]]@ranges@width,
                  ") does not match the chromosome length s_len in the blast_tbl (",
                  unique(species_specific_blast_tbl_plus_strand$s_len),
                  "). Thus, BLAST hit sequences cannot be extracted due to indexing shifts...",
                  call. = FALSE
                  
                )
              
              seqs_plus <- Biostrings::extractAt(
                imported_genome_i[chr_names[j]],
                at = IRanges::IRanges(
                  start = species_specific_blast_tbl_plus_strand$s_start,
                  end = species_specific_blast_tbl_plus_strand$s_end,
                  names = paste0(
                    species_specific_blast_tbl_plus_strand$species,
                    "_",
                    species_specific_blast_tbl_plus_strand$subject_id,
                    "_",
                    species_specific_blast_tbl_plus_strand$s_start,
                    "_",
                    species_specific_blast_tbl_plus_strand$s_end,
                    "_strand_plus"
                  )
                )
              )
              
              seqs_minus <- Biostrings::extractAt(
                Biostrings::reverse(imported_genome_i[chr_names[j]]),
                at = IRanges::IRanges(
                  start = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_start + 1,
                  end = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_end + 1,
                  names = paste0(
                    species_specific_blast_tbl_minus_strand$species,
                    "_",
                    species_specific_blast_tbl_minus_strand$subject_id,
                    "_",
                    species_specific_blast_tbl_minus_strand$s_start,
                    "_",
                    species_specific_blast_tbl_minus_strand$s_end,
                    "_strand_minus"
                  )
                )
              )
              
              Biostrings::writeXStringSet(
                seqs_plus@unlistData,
                filepath = ifelse(
                  is.null(path),
                  file.path(getwd(), file_name),
                  file.path(path, file_name)
                ),
                format = "fasta",
                append = TRUE
              )
              
              Biostrings::writeXStringSet(
                seqs_minus@unlistData,
                filepath = ifelse(
                  is.null(path),
                  file.path(getwd(), file_name),
                  file.path(path, file_name)
                ),
                format = "fasta",
                append = TRUE
              )
            }
          }
        }
      }
    }
    message("Sequence extraction process for all organisms finished without any problems.")
    return(seq_file_paths)
  }
