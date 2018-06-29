#' @title Perfrom BLAST Searches Against a Set of Genomes
#' @description This function takes a fasta file containing query sequences
#' as input and performs BLAST searches of these query sequences against
#' a set of reference genomes to retrieve corresponding hits in diverse genomes.
#' @param query path to input file in fasta format.
#' @param subject_genomes a vector containing file paths to the reference genomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).
#' @param blast_type specification of the BLAST type shall be used to perform BLAST searches between query and reference.
#' Available options are:
#' \itemize{
#' \item \code{blast_type = "blastn"} for  nucleotide to nucleotide BLAST searches (see \code{\link{blast_nucleotide_to_nucleotide}} for details)
#' \item \code{blast_type = "tblastn"} for protein to nucleotide BLAST searches (see \code{\link{blast_protein_to_nucleotide}} for details)
#' }
#' @param blast_output_path a path to a folder that will be created to store BLAST output tables for each individual query-genome search.
#' @param min_alig_length minimum alignment length that shall be retained in the result dataset. All hit alignments with smaller
#' hit alignment length will be removed automatically. 
#' @param \dots additional arguments passed to \code{\link{blast_nucleotide_to_nucleotide}} or \code{\link{blast_protein_to_nucleotide}}.
#' @author Hajk-Georg Drost
#' @details The \code{blast_genomes} function enables users to BLAST specific query sequences against a set of reference genomes
#' and retrieve the corresponding BLAST output.
#' @export


blast_genomes <-
  function(query,
           subject_genomes,
           blast_type = "blastn",
           blast_output_path = "blast_output",
           min_alig_length  = 50,
           ...) {
    
    if (!is.element(blast_type, c("blastn", "tblastn")))
      stop("Please specify a blast_type that is supported by this function, e.g. blast_type = 'blastn' or blast_type = 'tblastn'.", call. = FALSE)
    
    message("Start '",blast_type,"' search of query '", basename(query), "' against ", length(subject_genomes), " reference genomes ...")
    
    if (!file.exists(blast_output_path)) {
      message("BLAST results for each individual genome search will be stored at '", blast_output_path, "'.")
      dir.create(blast_output_path, recursive = TRUE)
    }
      
    
    res <- vector("list", length(subject_genomes))
    
    for (i in seq_len(length(subject_genomes))) {
      
      # remove *.fa appendix
      species_name <- unlist(stringr::str_split(basename(subject_genomes[i]),"[.]"))[1]
      
      message("Processing species: ", species_name," (",i, "/", length(subject_genomes),") ...")
      
      output_file <- paste0(species_name, "_genomes_blast.tbl")
      
      if (!file.exists(file.path(blast_output_path, output_file))) {
        
        if (blast_type == "blastn") {
          blast_output_tmp <- metablastr::blast_nucleotide_to_nucleotide(
            query = query,
            subject = subject_genomes[i],
            ...
          )
          if (!is.logical(blast_output_tmp)) {
            alig_length <- NULL
            blast_output_tmp <- dplyr::filter(blast_output_tmp, alig_length >= min_alig_length)
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, species = rep(species_name, nrow(blast_output_tmp)))
            
            message("Storing results for species ", species_name," in file ", file.path(blast_output_path, output_file), " ...")
            
            
            readr::write_excel_csv(blast_output_tmp,
                                   file.path(blast_output_path, output_file))
            
            res[i] <- blast_output_tmp
          }
        }
        
        if (blast_type == "tblastn") {
          blast_output_tmp <- metablastr::blast_protein_to_nucleotide(
            query = query,
            subject = subject_genomes[i],
            ...
          )
          if (!is.logical(blast_output_tmp)) {
            alig_length <- NULL
            blast_output_tmp <- dplyr::filter(blast_output_tmp, alig_length >= min_alig_length)
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, species = rep(species_name, nrow(blast_output_tmp)))
            
            message("Storing results for species ", species_name," in file ", file.path(blast_output_path, output_file), " ...")
            
            
            readr::write_excel_csv(blast_output_tmp,
                                   file.path(blast_output_path, output_file))
            
            res[i] <- blast_output_tmp
          }
        }
        
      } else {
        message("The blast output file '",file.path(blast_output_path, output_file),"' exists already and will be imported ...")
        res[i] <-
          readr::read_csv(file.path(blast_output_path, output_file))
      }
      
    }
    
    final_species_hit_tbl <- dplyr::bind_rows(res)
    message("The genome BLAST process has successfully been finished.")
    return(final_species_hit_tbl)
  }
