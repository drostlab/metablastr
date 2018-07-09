#' @title Extract sequences of BLAST hits in respective genomes and store it as \code{fasta} file(s)
#' @description In some cases, users may wish to extract the sequences of the respective blast hit
#' stored in a \code{blast_tbl}. This function enables to quickly extract such sequences and store
#' them in one general or multiple genome specific \code{fasta} file(s).
#' @param blast_tbl a BLAST table generated with \code{\link{blast_genomes}}.
#' @param subject_genomes a vector containing file paths to the reference genomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).
#' @param separated_by_genome a logical value indicating whether or not hit sequences from different genomes should be stored in the same
#' output \code{fasta} file \code{separated_by_genome = FALSE} (default) or in separate \code{fasta} files \code{separated_by_genome = TRUE}.
#' @param path a folder path in which corresponding \code{fasta} output files shall be stored.
#' @author Hajk-Georg Drost
#' @export
extract_hit_seqs_from_genomes <-
  function(blast_tbl,
           subject_genomes,
           separated_by_genome = FALSE,
           path = NULL) {
    
    if (any(!file.exists(subject_genomes)))
      stop(
        "At least one of the genome paths seems not to exist. Please check that all paths correspond to the correct location of the genome files.",
        call. = FALSE
      )
    
    message("Starting sequence extraction process ...")
    
    ifelse(separated_by_genome, seq_file_paths <- vector("character", 1), seq_file_paths <- vector("character", length(subject_genomes)))
    
    species <- s_strand <- NULL
    
    if (separated_by_genome) {
      for (i in seq_len(length(subject_genomes))) {
        message("Processing organism ", basename(subject_genomes[i]), " ...")
        imported_genome_i <- biomartr::read_genome(subject_genomes[i])
        
        # remove appendix *.fa from file name
        species_refined_name <- unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        species_specific_blast_tbl <-
          dplyr::filter(blast_tbl, species == species_refined_name)
        
        if (nrow(species_specific_blast_tbl) > 0) {
          strand <-
            ifelse(
              species_specific_blast_tbl$s_end - species_specific_blast_tbl$s_start + 1 >= 0,
              "+",
              "-"
            )
          
          species_specific_blast_tbl <-
            dplyr::mutate(species_specific_blast_tbl, s_strand = strand)
          
          # only retain chromosome names that are present in both: genome and BLAST table
          chr_names <- dplyr::intersect(names(imported_genome_i), names(table(species_specific_blast_tbl$subject_id)))
          
          if (length(chr_names) == 0)
            stop("It seems that the chromosome names used in the input 'blast_tbl' and in the input genome fasta file do not match. Please make sure that the chromosome names match in both cases.", call. = FALSE)

          # construct output file paths
          seq_file_paths[i] <- ifelse(
            is.null(path),
            file.path(getwd(), paste0(subject_genomes[i], "blast_tbl_sequences.fa")),
            file.path(path, paste0(subject_genomes[i], "blast_tbl_sequences.fa"))
          )
          
          # for each chromosome separately
          for (j in seq_len(length(chr_names))) {
            # divide by plus and minus strand
            species_specific_blast_tbl_plus_strand <-
              dplyr::filter(species_specific_blast_tbl, s_strand == "+")
            species_specific_blast_tbl_minus_strand <-
              dplyr::filter(species_specific_blast_tbl, s_strand == "-")
            
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
                start = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_end + 1,
                end = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_start + 1,
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
                file.path(getwd(), paste0(subject_genomes[i], "blast_tbl_sequences.fa")),
                file.path(path, paste0(subject_genomes[i], "blast_tbl_sequences.fa"))
              ),
              format = "fasta",
              append = TRUE
            )
            
            Biostrings::writeXStringSet(
              seqs_minus@unlistData,
              filepath = ifelse(
                is.null(path),
                file.path(getwd(), paste0(subject_genomes[i], "blast_tbl_sequences.fa")),
                file.path(path, paste0(subject_genomes[i], "blast_tbl_sequences.fa"))
              ),
              format = "fasta",
              append = TRUE
            )
          }
        } else {
          message("Organism ",
                  basename(subject_genomes[i]),
                  " didn't have any BLAST hits.")
        }
      }
    } else {
      seq_file_paths <- ifelse(
        is.null(path),
        file.path(getwd(), paste0("full_blast_tbl_sequences.fa")),
        file.path(path, paste0("full_blast_tbl_sequences.fa"))
      )
      
      for (i in seq_len(length(subject_genomes))) {
        message("Processing organism ", basename(subject_genomes[i]), " ...")
        imported_genome_i <- biomartr::read_genome(subject_genomes[i])
        
        # remove appendix *.fa from file name
        species_refined_name <- unlist(stringr::str_split(basename(subject_genomes[i]), "[.]"))[1]
        
        species_specific_blast_tbl <-
          dplyr::filter(blast_tbl, species == species_refined_name)
        
        if (nrow(species_specific_blast_tbl) > 0) {
          strand <-
            ifelse(
              species_specific_blast_tbl$s_end - species_specific_blast_tbl$s_start + 1 >= 0,
              "+",
              "-"
            )
          species_specific_blast_tbl <-
            dplyr::mutate(species_specific_blast_tbl, s_strand = strand)
          
          # only retain chromosome names that are present in both: genome and BLAST table
          chr_names <- dplyr::intersect(names(imported_genome_i), names(table(species_specific_blast_tbl$subject_id)))
          
          if (length(chr_names) == 0)
            stop("It seems that the chromosome names used in the input 'blast_tbl' and in the input genome fasta file do not match. Please make sure that the chromosome names match in both cases.", call. = FALSE)
          
          # for each chromosome separately
          for (j in seq_len(length(chr_names))) {
            # divide by plus and minus strand
            species_specific_blast_tbl_plus_strand <-
              dplyr::filter(species_specific_blast_tbl, s_strand == "+")
            species_specific_blast_tbl_minus_strand <-
              dplyr::filter(species_specific_blast_tbl, s_strand == "-")
            
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
                start = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_end + 1,
                end = imported_genome_i[chr_names[j]]@ranges@width - species_specific_blast_tbl_minus_strand$s_start + 1,
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
                file.path(getwd(), paste0("full_blast_tbl_sequences.fa")),
                file.path(path, paste0("full_blast_tbl_sequences.fa"))
              ),
              format = "fasta",
              append = TRUE
            )
            
            Biostrings::writeXStringSet(
              seqs_minus@unlistData,
              filepath = ifelse(
                is.null(path),
                file.path(getwd(), paste0("full_blast_tbl_sequences.fa")),
                file.path(path, paste0("full_blast_tbl_sequences.fa"))
              ),
              format = "fasta",
              append = TRUE
            )
          }
        } else {
          message("Organism ",
                  basename(subject_genomes[i]),
                  " didn't have any BLAST hits.")
        }
      }
    }
    message("Sequence extraction process finished properly.")
    return(seq_file_paths)
  }
