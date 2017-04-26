#' @title Perform Protein to Protein BLAST Searches
#' @description Run protein to protein BLAST of reference sequences
#' against a blast-able database.
#' @param query path to input file in fasta format.
#' @param subject path to subject file in fasta format or blast-able database.
#' @param is.subject.db logical specifying whether or not the \code{subject} file is a file in fasta format (\code{is.subject.db = FALSE}; default)
#' or a blast-able database that was formatted with \code{makeblastdb} (\code{is.subject.db = TRUE}).
#' @param task protein search task option. Options are:
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default).
#' \item \code{task = "blast-fast"} : Improved BLAST searches using longer words for protein seeding.
#' \item \code{task = "blastp-short"} : Optimized protein-protein comparisons for query sequences shorter than 30 residues.
#' }
#' @param import shall output of the protein BLAST search be directly imported via \code{\link{read_blast}}? Default is \code{import = FALSE}.
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 10}).
#' @param out.format a character string specifying the format of the file in which the BLAST results shall be stored.
#' Available options are:
#'  \itemize{
#'  \item \code{out.format = "pair"} : Pairwise
#'  \item \code{out.format = "qa.ident"} : Query-anchored showing identities
#'  \item \code{out.format = "qa.nonident"} : Query-anchored no identities
#'  \item \code{out.format = "fq.ident"} : Flat query-anchored showing identities
#'  \item \code{out.format = "fq.nonident"} : Flat query-anchored no identities
#'  \item \code{out.format = "xml"} : XML
#'  \item \code{out.format = "tab"} : Tabular separated file
#'  \item \code{out.format = "tab.comment"} : Tabular separated file with comment lines
#'  \item \code{out.format = "ASN.1.text"} : Seqalign (Text ASN.1)
#'  \item \code{out.format = "ASN.1.binary"} : Seqalign (Binary ASN.1)
#'  \item \code{out.format = "csv"} : Comma-separated values
#'  \item \code{out.format = "ASN.1"} : BLAST archive (ASN.1)
#'  \item \code{out.format = "json.seq.aln"} : Seqalign (JSON)
#'  \item \code{out.format = "json.blast.multi"} : Multiple-file BLAST JSON
#'  \item \code{out.format = "xml2.blast.multi"} : Multiple-file BLAST XML2
#'  \item \code{out.format = "json.blast.single"} : Single-file BLAST JSON
#'  \item \code{out.format = "xml2.blast.single"} : Single-file BLAST XML2
#'  \item \code{out.format = "report"} : Organism Report
#'  }
#' @param cores number of cores for parallel BLAST searches.
#' @param max.target.seqs maximum number of aligned sequences that shall be kept.
#' @param db.soft.mask shall low complexity regions be soft masked? Default is \code{db.soft.mask = FALSE}.
#' @param db.hard.mask shall low complexity regions be hard masked? Default is \code{db.hard.mask = FALSE}.
#' @param blast.path path to BLAST executables.
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' blast_protein(query   = system.file('seqs/qry.fa', package = 'metablastr'),
#'               subject = system.file('seqs/sbj.fa', package = 'metablastr'),
#'               import  = TRUE)
#' }
#' 
#' 
#' @export

blast_protein <- function(query, 
                          subject,
                          output.path = NULL,
                          is.subject.db = FALSE,
                          task = "blastp",
                          import = FALSE,
                          evalue   = 1E-3,
                          out.format = "tab", 
                          cores = 1,
                          max.target.seqs = 500,
                          db.soft.mask = FALSE,
                          db.hard.mask = FALSE,
                          blast.path = NULL) {
    
    if (!is_blast_installed())
        stop("Please install a valid version of BLAST.", call. = FALSE)
    
    # determine the number of cores on a multicore machine
    multi.cores <- parallel::detectCores()
    
    # in case one tries to use more cores than are available
    if (cores > multi.cores)
        stop("You chose more cores than are available on your machine.", call. = FALSE)
    
    # test if query file exists
    if (!file.exists(query))
        stop("Unfortunately, no query file has been found at ",query, call. = FALSE)
    
    # test if subject file or database exists
    if (!file.exists(subject))
        stop("Unfortunately, no subject file has been found at ", subject, call. = FALSE)
    
    if (!is.element(task, c("blastp", "blastp-fast", "blastp-short")))
        stop("Please choose a protein-protein comparison task that is supported by BLAST: task = 'blastp', task = 'blastp-fast', or task = 'blastp-short'.", call. = FALSE)
    
    blast_call <-
        paste0("blastp -query ", ws_wrap(query), " -db ", ws_wrap(subject))
    
    
    output_blast <-
        file.path(ifelse(is.null(output.path), ws_wrap(getwd()), ws_wrap(output.path)),
                  paste0(unlist(stringr::str_split(
                      basename(query), "[.]"
                  ))[1], ".blast_tbl"))
    
    output_read_blast <-
        file.path(ifelse(is.null(output.path), getwd(), output.path),
                  paste0(unlist(stringr::str_split(
                      basename(query), "[.]"
                  ))[1], ".blast_tbl"))
    
    params <-
        c(
            "query_id",
            "subject_id",
            "subject_taxonomy",
            "subject_kingdom",
            "perc_identity",
            "num_ident_matches",
            "alig_length",
            "mismatches",
            "gap_openings",
            "n_gaps",
            "pos_match",
            "ppos",
            "q_start",
            "q_end",
            "q_len",
            "qcov",
            "qcovhsp",
            "query_seq",
            "s_start",
            "s_end",
            "s_len",
            "subject_seq",
            "evalue",
            "bit_score",
            "score_raw"
        )
    
    # format subject into database
    if (!is.subject.db) {
        if (is.null(blast.path)) {
            system(
                paste0(
                    "makeblastdb -in ",
                    subject,
                    " -input_type fasta -dbtype prot -hash_index"
                )
            )
            
        } else {
            system(
                paste0(
                    "export PATH=",
                    blast.path,
                    "; makeblastdb -in ",
                    subject,
                    " -input_type fasta -dbtype prot -hash_index"
                )
            )
        }
    } 
    
    
    system(
        paste0(
            ifelse(is.null(blast.path), blast_call, paste0("export PATH=$PATH:", blast_call)),
            " -evalue ",
            evalue,
            " -max_target_seqs ",
            max.target.seqs,
            " -out ",
            output_blast ,
            " -num_threads ",
            cores,
            ifelse(db.soft.mask, " -db_soft_mask", ""),
            ifelse(db.hard.mask, " -db_hard_mask", ""),
            paste0( " -task ", task),
            paste0(" -outfmt '", outformat2num(out.format = out.format), " qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp sstart send slen evalue bitscore score'")
        )
    )
    
    
    if (import) {
        blast_tbl <- read_blast(file = output_read_blast, out.format = out.format)
        return(blast_tbl)
    }
    
}

