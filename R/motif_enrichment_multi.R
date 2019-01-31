#' @title 
#' @description 
#' @param 
#' @param 
#' @param 

motif_enrichment_multi <-
  function(blast_tbl,
           subject_genomes,
           test = "fisher",
           alternative = "two.sided",
           size,
           interval_width,
           motifs,
           max.mismatch = 0,
           min.mismatch = 0,
           ...) {
    species_names <- names(table(blast_tbl$species))
    res <- vector("list", length(species_names))
    species <- fisher_pval <- NULL

    for (i in seq_len(length(species_names))) {
      species_blast_tbl <-
        dplyr::filter(blast_tbl, species == species_names[i])
      species_subject_genome_path <-
        subject_genomes[stringr::str_detect(subject_genomes, species_names[i])]
      message(
        "Processing species ",
        species_names[i],
        " and subject file '",
        species_subject_genome_path,
        "' ..."
      )
      
      species_i_random_seqs <-
        metablastr::extract_random_seqs_from_genome(
          size = size,
          interval_width = interval_width,
          subject_genome = species_subject_genome_path,
          file_name = file.path(tempdir(), paste0(species_names[i], "_RandomSeqs.fa"))
        )
      
      extracted_seqs_path <- paste0(species_names[i], "_extracted_blast_hits.fa")
      
      metablastr::extract_hit_seqs_from_genomes(
        blast_tbl = species_blast_tbl,
        subject_genomes =  species_subject_genome_path,
        file_name = extracted_seqs_path,
        separated_by_genome = FALSE,
        path = tempdir()
      )
      
      species_motif_enrichment <- metablastr::motif_enrichment(
        real_seqs = file.path(tempdir(), extracted_seqs_path),
        random_seqs = species_i_random_seqs,
        test = test,
        alternative = alternative,
        motifs = motifs,
        max.mismatch = max.mismatch,
        min.mismatch = max.mismatch,
        ...
      )
      
      species_motif_enrichment <-
        dplyr::mutate(
          species_motif_enrichment,
          status = ifelse(fisher_pval <= 0.01, "significant", "not significant")
        )
      
      species_tbl <- tibble::tibble(species = rep(species_names[i], nrow(species_motif_enrichment)))
      
      species_motif_enrichment <-
        dplyr::bind_cols(species_tbl, species_motif_enrichment)
      
      res[i] <- list(species_motif_enrichment)
    }
    
    res <- dplyr::bind_rows(res)
    return(res)
  }


