#' @title Perform Protein to Protein BLAST Searches
#' @description Run protein to protein BLAST of reference sequences
#' against a blast-able database.
#' @param query path to input file in fasta format.
#' @param subject path to subject file in fasta format or blast-able database.
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
#'  },
#' @param cores
#' @param blast.path
#' @author Hajk-Georg Drost
#' @export
blast_p <- function(query, subject, out.format = "tab", cores = 1, blast.path = NULL) {
    if (!is_blast_installed())
        stop("Please install a valid version of BLAST.", call. = FALSE)
    
    # determine the number of cores on a multicore machine
    multi.cores <- parallel::detectCores()
    
    # in case one tries to use more cores than are available
    if (cores > multi.cores)
        stop("You chose more cores than are available on your machine.", call. = FALSE)
    
    if (is.null(blast.path)) {
                # use the default parameters when running blastp
        system(
            paste0(
                "blastp -db ",
                database,
                " -query ",
                input,
                " -evalue ",
                eval,
                " -max_target_seqs ",
                max.target.seqs,
                " -out ",
                output ,
                " -outfmt 6",
                " -num_threads ",
                cores
            )
        )
            } else {
                
                # use the default parameters when running blastp
                system(
                    paste0(
                        "export PATH=$PATH:",
                        path,
                        "; blastp -db ",
                        database,
                        " -query ",
                        input,
                        " -evalue ",
                        eval,
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
}

