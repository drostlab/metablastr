#' @title Retrieve all upstream promotor sequences from a genome
#' @description Given a genome assembly file and an corresponding annotation 
#' file users can retrieve all upstream promotor sequences of all genes from a genome. 
#' @param genome_file file path to the genome assembly file.
#' @param annotation_file file path to the annotation file of the genome assembly
#' in \code{gtf} or \code{gff} format.
#' @param promotor_length width of upstream promotors. This is -\code{promotor_width} bp from the transcription start site (TSS) of the gene.
#' @param annotation_format format of the annotation file. Options are:
#' \itemize{
#' \item \code{annotation_format = "gtf"}
#' \item \code{annotation_format = "gff"}
#' }
#' @param file_name file path to the output file storing the promotor sequences.
#' @param path a file path to an output folder storing the promotor sequences.  
#' @param update shall previously generated promotor sequences be overwritten when again generated for the same genome assembly?
#' @author Hajk-Georg Drost
#' @export
extract_promotor_seqs_from_genome <-
  function(annotation_file,
           genome_file,
           promotor_length = 500,
           annotation_format = "gtf",
           file_name = NULL,
           path = NULL,
           update = TRUE) {
    
    if (!file.exists(annotation_file))
      stop(
        "Please provide a valid path to the annotation file. The file '",
        annotation_file,
        "' seems not to exist.",
        call. = FALSE
      )
    
    if (!file.exists(genome_file))
      stop(
        "Please provide a valid path to the annotation file. The file '",
        annotation_file,
        "' seems not to exist.",
        call. = FALSE
      )
    
    if (promotor_length < 0)
      stop("A promotor length cannot be negative.", call. = FALSE)
    
    if (!is.element(annotation_format, c("gtf")))
      stop("Please specify an annotation_format that is supported by this function.", call. = FALSE)
    
    message(
      "Starting promotor sequence [length = ",
      promotor_length,
      "]",
      " extraction process of genome file ",
      genome_file,
      " and annotation file ",
      annotation_file,
      " ..."
    )
    
    message("Importing annotation file ...")
    
    
    if (annotation_format == "gtf") {
      Import_gtf <-
        tibble::as_tibble(rtracklayer::import(annotation_file))
    }
    
    type <- gene_biotype <- NULL
    
    message("Filtering annotation file for: type == 'gene', gene_biotype == 'protein_coding' ...")
    Import_gtf_filtered <-
      dplyr::filter(
        Import_gtf,
        type == "gene",
        gene_biotype == "protein_coding"
      )
    
    if (nrow(dplyr::filter(Import_gtf_filtered, (strand == "-") & (end - start + 1 <= 0))) > 0)
      stop("Your start and end coordinates for all minus strand genes are swapped (see: end - start + 1 <= 0).",
           " Please specify the minus strand as '-' in the strand column of your annotation file and use start and end coordinates analogous to the plus strand directtion, hence: end - start + 1 > 0.",call. = FALSE)
    
    message("Import genome file ...")
    Import_genome <- biomartr::read_genome(genome_file)
    
    if (nrow(Import_gtf_filtered) > 0) {
      Import_gtf_filtered_plus <- dplyr::filter(Import_gtf_filtered, (strand == "+") & (end - start + 1 > 0))
      Import_gtf_filtered_plus <- dplyr::mutate(Import_gtf_filtered_plus, strand_name = rep("plus", nrow(Import_gtf_filtered_plus)))
      
      Import_gtf_filtered_minus <- dplyr::filter(Import_gtf_filtered, (strand == "-") & (end - start + 1 > 0))
      Import_gtf_filtered_minus <- dplyr::mutate(Import_gtf_filtered_minus, strand_name = rep("minus", nrow(Import_gtf_filtered_minus)))
      
      Import_gtf_filtered <- dplyr::bind_rows(Import_gtf_filtered_plus, Import_gtf_filtered_minus)
      
      # only retain chromosome names that are present in both: genome and BLAST table
      chr_names <-
        dplyr::intersect(unique(Import_genome@ranges@NAMES), names(table(Import_gtf_filtered$seqnames)))
      chr_names <- as.character(chr_names)
      
      if (length(chr_names) == 0) {
        Import_genome@ranges@NAMES <-
          as.character(sapply(Import_genome@ranges@NAMES, function(x)
            unlist(stringr::str_split(x, " "))[1]))
        
        chr_names <-
          dplyr::intersect(unique(Import_genome@ranges@NAMES), names(table(Import_gtf_filtered$seqnames)))
        
        if (length(chr_names) == 0) {
          stop(
            "It seems that the chromosome names used in the input  annotation file and in the input genome fasta file do not match. Please make sure that the chromosome names match in both cases.",
            call. = FALSE
          )
        }
      }
      
      # construct output file paths
      promotor_seq_file_output <- ifelse(is.null(path),
                                         file.path(getwd(),
                                                   paste0(
                                                     unlist(stringr::str_split(basename(genome_file), "[.]"))[1],
                                                     paste0("_promotor_seqs_genes_length_", promotor_length)
                                                   )),
                                         file.path(path,
                                                   paste0(
                                                     unlist(stringr::str_split(basename(genome_file), "[.]"))[1],
                                                     paste0("_promotor_seqs_genes_length_", promotor_length)
                                                   )))
      
      if (update & file.exists(promotor_seq_file_output))
        fs::file_delete(promotor_seq_file_output)
      
      # for each chromosome separately
      for (j in seq_len(length(chr_names))) {
        message("  -> Extracting promotor sequences from ",
                chr_names[j],
                " ...")
        
        seqnames <- strand_name <- NULL
        # divide by plus and minus strand
        ## plus strand
        Import_gtf_filtered_plus_strand <-
          dplyr::filter(Import_gtf_filtered,
                        seqnames == chr_names[j],
                        strand_name == "plus")
        
        start <- end <- NULL
        Import_gtf_filtered_plus_strand_pissible <- dplyr::filter(Import_gtf_filtered_plus_strand, start - promotor_length > 0)
        
        if (nrow(Import_gtf_filtered_plus_strand) - nrow(Import_gtf_filtered_plus_strand_pissible) > 0)
          message("There were ",nrow(Import_gtf_filtered_plus_strand) - nrow(Import_gtf_filtered_plus_strand_pissible), " gene(s) on the plus strand for which promotor sequences couldn't be extracted because the promotor regions were beyond the borders of ", chr_names[j] , " ...")
        
        ## minus strand
        Import_gtf_filtered_minus_strand <-
          dplyr::filter(Import_gtf_filtered,
                        seqnames == chr_names[j],
                        strand_name == "minus")
        
        Import_gtf_filtered_minus_strand_pissible <- dplyr::filter(Import_gtf_filtered_minus_strand, end - promotor_length > 0)
        
        if (nrow(Import_gtf_filtered_minus_strand) - nrow(Import_gtf_filtered_minus_strand_pissible) > 0)
          message("There were ", nrow(Import_gtf_filtered_minus_strand) - nrow(Import_gtf_filtered_minus_strand_pissible), " genes on the minus strand for which promotor sequences couldn't be extracted because the promotor regions were beyond the borders of ", chr_names[j] , " ...")
        
        if (nrow(Import_gtf_filtered_plus_strand_pissible) > 0) {
          seqs_plus <- Biostrings::extractAt(
            Import_genome[chr_names[j]],
            at = IRanges::IRanges(
              start = Import_gtf_filtered_plus_strand_pissible$start - promotor_length,
              end = Import_gtf_filtered_plus_strand_pissible$start - 1,
              names = paste0(
                Import_gtf_filtered_plus_strand_pissible$gene_id,
                "_",
                Import_gtf_filtered_plus_strand_pissible$seqnames,
                "_",
                Import_gtf_filtered_plus_strand_pissible$start - promotor_length,
                "_",
                Import_gtf_filtered_plus_strand_pissible$start - 1,
                "_strand_plus"
              )
            )
          )
          
          # store sequences in fasta files
          Biostrings::writeXStringSet(
            seqs_plus@unlistData,
            filepath = promotor_seq_file_output,
            format = "fasta",
            append = TRUE
          )
        } else {
          warning("The plus strand of ",chr_names[j], " did not contain any genes and thus promotor extraction was omitted here.", call. = FALSE)
        }
        
        if (nrow(Import_gtf_filtered_minus_strand_pissible) > 0) {
          seqs_minus <-
            Biostrings::reverseComplement(Biostrings::extractAt(
              Import_genome[chr_names[j]],
              at = IRanges::IRanges(
                start = Import_gtf_filtered_minus_strand_possible$end + 1,
                end = Import_gtf_filtered_minus_strand_possible$end + promotor_length,
                names = paste0(
                  Import_gtf_filtered_minus_strand_possible$gene_id,
                  "_",
                  Import_gtf_filtered_minus_strand_possible$seqnames,
                  "_",
                  Import_gtf_filtered_minus_strand_possible$end - promotor_length,
                  "_",
                  Import_gtf_filtered_minus_strand_possible$end -  1,
                  "_strand_minus"
                )
              )
            )[[1]])
          
          Biostrings::writeXStringSet(seqs_minus,
                                      filepath = promotor_seq_file_output,
                                      format = "fasta",
                                      append = TRUE)

        } else {
          warning("The minus strand of ",chr_names[j], " did not contain any genes and thus promotor extraction was omitted here.", call. = FALSE)
        }
      }
    } else {
      stop("After filtering for genes the annotation file was empty. Please check what might have gone wrong!", call. = FALSE)
    }
    
    fs::file_move(promotor_seq_file_output, paste0(promotor_seq_file_output, ".fa"))
    
    message("Sequence extraction process for all organisms finished without any problems and output file was stored at ",paste0(promotor_seq_file_output, ".fa"),".")
    return(paste0(promotor_seq_file_output, ".fa"))
  }