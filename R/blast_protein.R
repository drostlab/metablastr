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
#' @param blast.path path to BLAST executables.
#' @author Hajk-Georg Drost
#' @export

blast_protein <- function(query, 
                          subject,
                          is.subject.db = FALSE,
                          task = "blastp",
                          import = FALSE,
                          evalue   = 10,
                          out.format = "tab", 
                          cores = 1, 
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
    
    ifelse(!is.subject.db, blast_call <- paste0(task, " -query ", query, " -subject ", subject),
                           blast_call <- paste0(task, " -query ", query, " -db ", subject))
    
    system(
        paste0(
            ifelse(is.null(blast.path), blast_call, paste0("export PATH=$PATH:", blast_call)),
            " -evalue ",
            evalue,
            " -max_target_seqs ",
            max.target.seqs,
            " -out ",
            output ,
            " -outfmt 6",
            " -num_threads ",
            cores
        )
    )
    
}

