#' @title Count the number of motifs in a set of non-random versus randomly sampled genomic sequences for multiple species
#' @description Compare the number of motifs in a set of non-random versus random sequences within a set of subject genomes. 
#' The resulting values can then be used to test the enrichment of certain motifs in real sequences compared to random sequences. 
#' Each enrichment analysis is performed for a set of different species or genomes.
#' @param blast_tbl a blast_table.
#' @param subject_genomes a character vector storing the file paths to the subject genomes that shall be used as subject references.
#' @param size total number of sequences that shall be sampled per subject genome.
#' @param interval_width length of the sequence in which motifs shall be detected.
#' @param motifs a character vector storing (case sensitive) motif sequences for which abundance in the sampled sequences shall be assessed.
#' @param max.mismatch maximum number of mismatches that are allowed between the sequence motif and the matching region in the sampled sequence.
#' @param min.mismatch minimum number of mismatches that are allowed between the sequence motif and the matching region in the sampled sequence.
#' @param ... additional arguments passed to \code{\link{motif_compare}}.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{motif_count}}, \code{\link{motif_compare_multi_promotor_seqs}}, \code{\link{motif_compare}}, 
#' \code{\link{motif_enrichment}}, \code{\link{motif_enrichment_multi}}
#' @export 

motif_compare_multi <-
  function(blast_tbl,
           subject_genomes,
           size,
           interval_width,
           motifs,
           max.mismatch = 0,
           min.mismatch = 0,
           ...) {
    species_names <- names(table(blast_tbl$species))
    res <- vector("list", length(species_names))
    species <- NULL 
    
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
      
      species_motif_compare <- metablastr::motif_compare(
        real_seqs = file.path(tempdir(), extracted_seqs_path),
        random_seqs = species_i_random_seqs,
        motifs = motifs,
        max.mismatch = max.mismatch,
        min.mismatch = max.mismatch,
        ...
      )
      
      species_motif_compare <-
        dplyr::mutate(species_motif_compare, species = rep(species_names[i], nrow(species_motif_compare)))
      species_motif_compare <-
        dplyr::select(species_motif_compare, species, 1:(ncol(species_motif_compare) - 1))
      res[i] <- list(species_motif_compare)
    }
    
    res <- dplyr::bind_rows(res)
    return(res)
  }
