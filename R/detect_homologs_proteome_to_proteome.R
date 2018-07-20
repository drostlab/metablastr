#' @title Perform Proteome to Proteome BLAST Searches against a set of Proteomes
#' @description Run proteome to proteome BLAST searches to detect homologous protein sequences
#' in a set of subject proteomes.
#' @param query path to input file in fasta format.
#' @param subject_proteomes a character vector containing paths to subject files in fasta format.
#' @param task protein search task option. Options are:
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default).
#' \item \code{task = "blast-fast"} : Improved BLAST searches using longer words for protein seeding.
#' \item \code{task = "blastp-short"} : Optimized protein-protein comparisons for query sequences shorter than 30 residues.
#' }
#' @param blast_output_path a path to a folder that will be created to store BLAST output tables for each individual query-proteome search.
#' @param min_alig_length minimum alignment length that shall be retained in the result dataset. All hit alignments with smaller hit alignment length will be removed automatically.
#' @param evalue Expectation value (E) threshold for saving hits (default: evalue = 1E-5).
#' @param max.target.seqs maximum number of aligned sequences that shall be kept. Default is \code{max.target.seqs = 500}.
#' @param cores number of cores for parallel BLAST searches.
#' @param update a logical value indicating whether or not pre-computed BLAST tables should be removed and re-computed (\code{update = TRUE}) or imported from existing file (\code{update = FALSE}) (Default).
#' @param \dots additional arguments passed to \code{\link{blast_protein_to_protein}}.
#' @author Hajk-Georg Drost
#' @export

detect_homologs_proteome_to_proteome <-
  function(query,
           subject_proteomes,
           task = "blastp",
           blast_output_path = "blast_output",
           min_alig_length  = 20,
           evalue = 1E-5,
           max.target.seqs = 5000,
           cores = 1,
           update = FALSE,
           ...) {
    if (!all(file.exists(subject_proteomes)))
      stop(
        "It seems that at least one proteome file path does not exist. Please check that all subject_proteomes comprise file paths leading to a .fasta file containing protein sequences.",
        call. = FALSE
      )
    
    message(
      "Start '",
      task,
      "' search of query '",
      basename(query),
      "' against ",
      length(subject_proteomes),
      " reference proteomes using evalue = ",
      evalue,
      " and max.target.seqs = ",
      max.target.seqs,
      " ..."
    )
    
    if (update) {
      file.remove(blast_output_path, recursive = TRUE, overwrite = TRUE)
    }
    
    if (!file.exists(blast_output_path)) {
      message(
        "BLAST results for each individual proteome search will be stored at '",
        blast_output_path,
        "'."
      )
      dir.create(blast_output_path, recursive = TRUE)
    }
    
    res <- vector("list", length(subject_proteomes))
    
    for (i in seq_len(length(subject_proteomes))) {
      # remove *.fa appendix
      species_name <-
        unlist(stringr::str_split(basename(subject_proteomes[i]), "[.]"))[1]
      
      message("Processing species: ",
              species_name,
              " (",
              i,
              "/",
              length(subject_proteomes),
              ") ...")
      
      output_file <- paste0(species_name, "_proteome_blast.tbl")
      
      if (!file.exists(file.path(blast_output_path, output_file))) {
        blast_query_vs_subject_i <- blast_protein_to_protein(
          query = query,
          subject = subject_proteomes[i],
          task = task,
          evalue = evalue,
          max.target.seqs = max.target.seqs,
          cores = cores,
          ...
        )
        
        if (!is.logical(blast_query_vs_subject_i)) {
          alig_length <- q_len <- NULL
          blast_query_vs_subject_i <-
            dplyr::filter(blast_query_vs_subject_i,
                          alig_length >= min_alig_length)
          blast_query_vs_subject_i <-
            dplyr::mutate(blast_query_vs_subject_i,
                          species = rep(species_name, nrow(blast_query_vs_subject_i)))
          blast_query_vs_subject_i <-
            dplyr::mutate(blast_query_vs_subject_i, scope = 1 - (abs(q_len - alig_length) / q_len))
          
          message(
            "Storing results for species ",
            species_name,
            " in file ",
            file.path(blast_output_path, output_file),
            " ..."
          )
          
          readr::write_excel_csv(blast_query_vs_subject_i,
                                 file.path(blast_output_path, output_file))
          
          res[i] <- list(blast_query_vs_subject_i)
        }
        
      } else {
        message(
          "The blast output file '",
          file.path(blast_output_path, output_file),
          "' exists already and will be imported ..."
        )
        suppressMessages(res[i] <-
                           list(readr::read_tsv(
                             file.path(blast_output_path, output_file)
                           )))
      }
    }
    
    final_species_hit_tbl <- dplyr::bind_rows(res)
    message("The proteome BLAST process has successfully been finished.")
    return(final_species_hit_tbl)
  }
