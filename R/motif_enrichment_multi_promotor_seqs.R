#' @title Statistical assessment of motif enrichments in a set of non-random versus randomly sampled gene promotor sequences for multiple species
#' @description Compare the number of motifs in a set of non-random versus randomly sampled gene promotor sequences within a set of subject genomes. 
#' The resulting values are then used to statistically assess the enrichment of certain motifs in real sequences compared to randomly sampled gene promotor sequences. 
#' @param blast_tbl a blast_table.
#' @param subject_genomes a character vector storing the file paths to the subject genomes that shall be used as subject references.
#' @param annotation_files a character vector storing the file paths to the subject annotation files in \code{.gff} format that match the subject genomes.
#' @param annotation_format the annotation format. Options are:
#' \itemize{
#' \item \code{annotation_format = "gff"}
#' }
#' @param test \itemize{
#' \item \code{test = "fisher"}: Fisher's Exact Test for Count Data (see \code{link[stats]{fisher.test}} for details).
#' }
#' @param alternative indicates the alternative hypothesis and must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just the initial letter. Only used in the 2 by 2 case.
#' @param interval_width total number of sequences that shall be sampled per subject genome.
#' @param motifs a character vector storing (case sensitive) motif sequences for which abundance in the sampled sequences shall be assessed.
#' @param max.mismatch maximum number of mismatches that are allowed between the sequence motif and the matching region in the sampled sequence.
#' @param min.mismatch minimum number of mismatches that are allowed between the sequence motif and the matching region in the sampled sequence.
#' @param ... additional arguments passed to \code{\link{motif_compare}}.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{motif_count}}, \code{\link{motif_compare_multi_promotor_seqs}}, \code{\link{motif_compare}}, 
#' \code{\link{motif_enrichment}}

motif_enrichment_multi_promotor_seqs <-
  function(blast_tbl,
           subject_genomes,
           annotation_files,
           annotation_format = "gff",
           test = "fisher",
           alternative = "two.sided",
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
        extract_upstream_promotor_seqs(
          organism = species_names[i],
          annotation_file = annotation_files[i],
          annotation_format = annotation_format,
          promotor_width = interval_width,
          genome_file = species_subject_genome_path,
          file_name = file.path(tempdir(), paste0(species_names[i], "_PromotorSeqs.fa")))
      
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


