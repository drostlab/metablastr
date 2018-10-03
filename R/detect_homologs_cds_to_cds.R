#' @title Perform CDS to CDS BLAST Searches against a set of CDS files
#' @description Run cds to cds BLAST searches to detect homologous cds sequences
#' in a set of subject cds files.
#' @param query path to input file in fasta format.
#' @param subject_proteomes a character vector containing paths to subject files in fasta format.
#' @param task nucleotide search task option. Options are:
#' \itemize{
#' \item \code{task = "blastn"} : Standard nucleotide-nucleotide comparisons (default) - Traditional BLASTN requiring an exact match of 11.
#' \item \code{task = "blastn-short"} : Optimized nucleotide-nucleotide comparisons for query sequences shorter than 50 nucleotides.
#' \item \code{task = "dc-megablast"} : Discontiguous megablast used to find more distant (e.g., interspecies) sequences.
#' \item \code{task = "megablast"} : Traditional megablast used to find very similar (e.g., intraspecies or closely related species) sequences.
#' \item \code{task = "rmblastn"}
#' }
#' @param blast_output_path a path to a folder that will be created to store BLAST output tables for each individual query-cds search.
#' @param min_alig_length minimum alignment length that shall be retained in the result dataset. All hit alignments with smaller hit alignment length will be removed automatically.
#' @param evalue Expectation value (E) threshold for saving hits (default: evalue = 1E-5).
#' @param max.target.seqs maximum number of aligned sequences that shall be kept. Default is \code{max.target.seqs = 500}.
#' @param cores number of cores for parallel BLAST searches.
#' @param update a logical value indicating whether or not pre-computed BLAST tables should be removed and re-computed (\code{update = TRUE}) or imported from existing file (\code{update = FALSE}) (Default).
#' @param \dots additional arguments passed to \code{\link{blast_nucleotide_to_nucleotide}}.
#' @author Hajk-Georg Drost
#' @export

detect_homologs_cds_to_cds <-
  function(query,
           subject_cds,
           task = "blastn",
           blast_output_path = "blast_output",
           min_alig_length  = 60,
           evalue = 1E-5,
           max.target.seqs = 5000,
           cores = 1,
           update = FALSE,
           ...) {
    if (!all(file.exists(subject_proteomes)))
      stop(
        "It seems that at least one cds file path does not exist. Please check that all subject_cds comprise file paths leading to a .fasta file containing cds sequences.",
        call. = FALSE
      )
    
    message(
      "Start '",
      task,
      "' search of query '",
      basename(query),
      "' against ",
      length(subject_proteomes),
      " reference cds using evalue = ",
      evalue,
      " and max.target.seqs = ",
      max.target.seqs,
      " ..."
    )
    
    if (update) {
      fs::dir_delete(blast_output_path)
    }
    
    if (!file.exists(blast_output_path)) {
      message(
        "BLAST results for each individual proteome search will be stored at '",
        blast_output_path,
        "'."
      )
      dir.create(blast_output_path, recursive = TRUE)
    }
    
    res <- vector("list", length(subject_cds))
    
    for (i in seq_len(length(subject_cds))) {
      # remove *.fa appendix
      species_name <-
        unlist(stringr::str_split(basename(subject_cds[i]), "[.]"))[1]
      
      message("Processing species: ",
              species_name,
              " (",
              i,
              "/",
              length(subject_cds),
              ") ...")
      
      output_file <- paste0(species_name, "_proteome_blast.tsv")
      
      if (!file.exists(file.path(blast_output_path, output_file))) {
        blast_query_vs_subject_i <- blast_nucleotide_to_nucleotide(
          query = query,
          subject = subject_cds[i],
          task = task,
          evalue = evalue,
          max.target.seqs = max.target.seqs,
          cores = cores,
          ...
        )
        
        if (file.exists(paste0(subject_cds[i],".nhd")))
          file.remove(paste0(subject_cds[i],".nhd"))  
        if (file.exists(paste0(subject_cds[i],".nhi")))
          file.remove(paste0(subject_cds[i],".nhi"))  
        if (file.exists(paste0(subject_cds[i],".nhr")))
          file.remove(paste0(subject_cds[i],".nhr"))
        if (file.exists(paste0(subject_cds[i],".nin")))
          file.remove(paste0(subject_cds[i],".nin")) 
        if (file.exists(paste0(subject_cds[i],".nog")))
          file.remove(paste0(subject_cds[i],".nog"))
        if (file.exists(paste0(subject_cds[i],".nog")))
          file.remove(paste0(subject_cds[i],".nog"))
        if (file.exists(paste0(subject_cds[i],".nsd")))
          file.remove(paste0(subject_cds[i],".nsd"))
        if (file.exists(paste0(subject_cds[i],".nsi")))
          file.remove(paste0(subject_cds[i],".nsi"))
        if (file.exists(paste0(subject_cds[i],".nsq")))
          file.remove(paste0(subject_cds[i],".nsq"))
        
        
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
          
          readr::write_tsv(blast_query_vs_subject_i,
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
    message("The cds BLAST process has successfully been finished.")
    return(final_species_hit_tbl)
  }
