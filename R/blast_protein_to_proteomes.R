#' @title Perfrom BLAST Searches Against a Set of Proteomes
#' @description This function takes a fasta file containing query sequences
#' as input and performs BLAST searches of these query sequences against
#' a set of reference proteomes to retrieve corresponding hits in diverse proteomes
#' @param query path to input file in fasta format.
#' @param subject_proteomes a vector containing file paths to the reference proteomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).
#' @param blast_type specification of the BLAST type shall be used to perform BLAST searches between query and reference.
#' Available options are:
#' \itemize{
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default); (see \code{\link{blast_protein_to_protein}} for details).
#' \item \code{task = "blast-fast"} : Improved BLAST searches using longer words for protein seeding.
#' \item \code{task = "blastp-short"} : Optimized protein-protein comparisons for query sequences shorter than 30 residues.
#' }
#' @param blast_output_path a path to a folder that will be created to store BLAST output tables for each individual query-proteome search.
#' @param min_alig_length minimum alignment length that shall be retained in the result dataset. All hit alignments with smaller
#' hit alignment length will be removed automatically.
#' @param evalue minimum expectation value (E) threshold for retaining hits (default: evalue = 0.00001).
#' @param max.target.seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param update a logical value indicating whether or not pre-computed BLAST tables should be removed and re-computed (\code{update = TRUE})
#' or imported from existing file (\code{update = FALSE}) (Default). 
#' @param \dots additional arguments passed to \code{\link{blast_protein_to_protein}}.
#' @author Hajk-Georg Drost
#' @details The \code{blast_protein_to_proteomes} function enables users to BLAST specific query sequences against a set of reference proteomes
#' and retrieve the corresponding BLAST output.
#' @export

blast_protein_to_proteomes <-
  function(query,
           subject_proteomes,
           blast_type = "blastp",
           blast_output_path = "blast_output",
           min_alig_length  = 30,
           evalue = 1E-5,
           max.target.seqs = 5000,
           update = FALSE,
           ...) {
    
    if (!is.element(blast_type, c("blastp", "blast-fast", "blastp-short")))
      stop("Please specify a blast_type that is supported by this function, e.g. blast_type = 'blastp',  blast_type = 'blast-fast',  or blast_type = 'blastp-short'.", call. = FALSE)
    
    message("Start '",blast_type,"' search of query '", basename(query), "' against ", length(subject_proteomes), " reference proteomes using evalue = ",evalue," and max.target.seqs = ",max.target.seqs," ...")
    
    if (update) {
      fs::dir_delete(blast_output_path)
    }
    
    if (!file.exists(blast_output_path)) {
      message("BLAST results for each individual proteome search will be stored at '", blast_output_path, "'.")
      dir.create(blast_output_path, recursive = TRUE)
    }
    
    res <- vector("list", length(subject_proteomes))
    
    for (i in seq_len(length(subject_proteomes))) {
      
      # remove *.fa appendix
      species_name <- unlist(stringr::str_split(basename(subject_proteomes[i]),"[.]"))[1]
      
      message("Processing species: ", species_name," (",i, "/", length(subject_proteomes),") ...")
      
      output_file <- paste0(species_name, "_proteomes_blast.tbl")
      
      if (!file.exists(file.path(blast_output_path, output_file))) {
        
        if (blast_type == "blastp") {
          blast_output_tmp <- metablastr::blast_protein_to_protein(
            query = query,
            subject = subject_proteomes[i],
            evalue = evalue,
            max.target.seqs = max.target.seqs,
            ...
          )
          
          # remove all makeblastdb files                    
          if (file.exists(paste0(subject_proteomes[i],".phd")))
            file.remove(paste0(subject_proteomes[i],".phd"))  
          if (file.exists(paste0(subject_proteomes[i],".phi")))
            file.remove(paste0(subject_proteomes[i],".phi"))  
          if (file.exists(paste0(subject_proteomes[i],".phr")))
            file.remove(paste0(subject_proteomes[i],".phr"))
          if (file.exists(paste0(subject_proteomes[i],".pin")))
            file.remove(paste0(subject_proteomes[i],".pin")) 
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".psd")))
            file.remove(paste0(subject_proteomes[i],".psd"))
          if (file.exists(paste0(subject_proteomes[i],".psi")))
            file.remove(paste0(subject_proteomes[i],".psi"))
          if (file.exists(paste0(subject_proteomes[i],".psq")))
            file.remove(paste0(subject_proteomes[i],".psq"))
          
          if (!is.logical(blast_output_tmp)) {
            alig_length <- q_len <- NULL
            blast_output_tmp <- dplyr::filter(blast_output_tmp, alig_length >= min_alig_length)
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, species = rep(species_name, nrow(blast_output_tmp)))
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, scope = 1 - (abs(q_len - alig_length) / q_len))
            
            message("Storing results for species ", species_name," in file ", file.path(blast_output_path, output_file), " ...")
            
            readr::write_excel_csv(blast_output_tmp,
                                   file.path(blast_output_path, output_file))
            
            res[i] <- list(blast_output_tmp)
          }
        }
        
        if (blast_type == "blast-fast") {
          blast_output_tmp <- metablastr::blast_protein_to_protein(
            query = query,
            subject = subject_proteomes[i],
            evalue = evalue,
            max.target.seqs = max.target.seqs,
            task = "blast-fast",
            ...
          )
          # remove all makeblastdb files                    
          if (file.exists(paste0(subject_proteomes[i],".phd")))
            file.remove(paste0(subject_proteomes[i],".phd"))  
          if (file.exists(paste0(subject_proteomes[i],".phi")))
            file.remove(paste0(subject_proteomes[i],".phi"))  
          if (file.exists(paste0(subject_proteomes[i],".phr")))
            file.remove(paste0(subject_proteomes[i],".phr"))
          if (file.exists(paste0(subject_proteomes[i],".pin")))
            file.remove(paste0(subject_proteomes[i],".pin")) 
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".psd")))
            file.remove(paste0(subject_proteomes[i],".psd"))
          if (file.exists(paste0(subject_proteomes[i],".psi")))
            file.remove(paste0(subject_proteomes[i],".psi"))
          if (file.exists(paste0(subject_proteomes[i],".psq")))
            file.remove(paste0(subject_proteomes[i],".psq"))
          
          if (!is.logical(blast_output_tmp)) {
            alig_length <- q_len <- NULL
            blast_output_tmp <- dplyr::filter(blast_output_tmp, alig_length >= min_alig_length)
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, species = rep(species_name, nrow(blast_output_tmp)))
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, scope = 1 - (abs(q_len - alig_length) / q_len))
            
            message("Storing results for species ", species_name," in file ", file.path(blast_output_path, output_file), " ...")
            
            readr::write_tsv(blast_output_tmp,
                             file.path(blast_output_path, output_file))
            
            res[i] <- list(blast_output_tmp)
          }
        }
        
        if (blast_type == "blastp-short") {
          blast_output_tmp <- metablastr::blast_protein_to_protein(
            query = query,
            subject = subject_proteomes[i],
            evalue = evalue,
            max.target.seqs = max.target.seqs,
            task = "blastp-short",
            ...
          )
          # remove all makeblastdb files                    
          if (file.exists(paste0(subject_proteomes[i],".phd")))
            file.remove(paste0(subject_proteomes[i],".phd"))  
          if (file.exists(paste0(subject_proteomes[i],".phi")))
            file.remove(paste0(subject_proteomes[i],".phi"))  
          if (file.exists(paste0(subject_proteomes[i],".phr")))
            file.remove(paste0(subject_proteomes[i],".phr"))
          if (file.exists(paste0(subject_proteomes[i],".pin")))
            file.remove(paste0(subject_proteomes[i],".pin")) 
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".pog")))
            file.remove(paste0(subject_proteomes[i],".pog"))
          if (file.exists(paste0(subject_proteomes[i],".psd")))
            file.remove(paste0(subject_proteomes[i],".psd"))
          if (file.exists(paste0(subject_proteomes[i],".psi")))
            file.remove(paste0(subject_proteomes[i],".psi"))
          if (file.exists(paste0(subject_proteomes[i],".psq")))
            file.remove(paste0(subject_proteomes[i],".psq"))
          
          if (!is.logical(blast_output_tmp)) {
            alig_length <- q_len <- NULL
            blast_output_tmp <- dplyr::filter(blast_output_tmp, alig_length >= min_alig_length)
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, species = rep(species_name, nrow(blast_output_tmp)))
            blast_output_tmp <- dplyr::mutate(blast_output_tmp, scope = 1 - (abs(q_len - alig_length) / q_len))
            
            message("Storing results for species ", species_name," in file ", file.path(blast_output_path, output_file), " ...")
            
            readr::write_tsv(blast_output_tmp,
                             file.path(blast_output_path, output_file))
            
            res[i] <- list(blast_output_tmp)
          }
        }
      } else {
        message("The blast output file '",file.path(blast_output_path, output_file),"' exists already and will be imported ...")
        suppressMessages(res[i] <-
                           list(readr::read_tsv(file.path(blast_output_path, output_file))))
      }
    }
    
    final_species_hit_tbl <- dplyr::bind_rows(res)
    message("The proteome BLAST process has successfully been finished.")
    return(final_species_hit_tbl)
  }
