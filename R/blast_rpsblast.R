#' @title Perform Reverse PSI-BLAST searches (rpsblast)
#' @description Run rpsblast (Reverse PSI-BLAST) searches a query sequence against a database of profiles, or score matrices, producing BLAST-like output.
#' @param query path to input file in fasta format.
#' @param db path to rpsblast-able database.
#' @param db.alias alias for database files 
#' @param output.path path to folder at which BLAST output table shall be stored. 
#' @param db.import shall the BLAST output be stored in a PostgresSQL database and shall a connection be established to this database? Default is \code{db.import = FALSE}.
#' In case users wish to to only generate a BLAST output file without importing it to the current R session they can specify \code{db.import = NULL}.
#' @param postgres.user when \code{db.import = TRUE} and \code{out.format = "postgres"} is selected, the BLAST output is imported and stored in a 
#' PostgresSQL database. In that case, users need to have PostgresSQL installed and initialized on their system. 
#' Please consult the Installation Vignette for details. 
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 0.001}).
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
#'  \item \code{out.format = "SAM"} : Sequence Alignment/Map (SAM)
#'  \item \code{out.format = "report"} : Organism Report
#'  }
#' @param cores number of cores for parallel rpsblast searches.
#' @param blast.path path to BLAST executables.
#' @author Anna Gogleva
#' @seealso \code{\link{blast_protein_to_protein}}, \code{\link{blast_nucleotide_to_protein}}
#' @export

blast_rpsblast <- function(query,
                           db,
                           db.alias,
                           output.path = NULL,
                           db.import = FALSE,
                           postgress.user = NULL,
                           evalue = 1E-3,
                           out.format = 'csv',
                           cores = 1,
                           blast.path = NULL) {
    
    if (!is_blast_installed())
        stop("Please install a valid version of rpsblast. See Installation Vignette for details.", call. = FALSE)
    
    if (db.import) {
        if (!is.element(out.format, c("xml", "tab", "csv")))
            stop("Only output formats: 'xml', 'tab', or 'csv' can be imported.", call. = FALSE)
    }
    
    # determine the number of cores on a multicore machine
    multi.cores <- parallel::detectCores()
    
    # in case one tries to use more cores than are available
    if (cores > multi.cores)
        stop("You chose more cores than are available on your machine.", call. = FALSE)
    
    # test if query file exists
    if (!file.exists(query))
        stop("Unfortunately, no query file has been found at ", query, call. = FALSE)
    
    # test if all required database files exist
    db_postfix <- c('aux', 'freq', 
                    'loo', 'phr',
                    'pin', 'psi',
                    'psq', 'rps')

    full_path <- paste(db, db.alias, sep = '/')
    db_files <- sapply(db_postfix,
                       function(x) paste(full_path, x, sep = '.'))
    db_status <- sapply(db_files, file.exists)
    
    if (!(all(db_stats))) {
        stop('Unfortunately, one ore more required database files are missing', names(db_status[db_status == FALSE]) call. = FALSE)
    }
    
    message("Starting 'rpsblast","' with  query: ", query, " against database: ", db.alias," using ", cores, " core(s) ...")
    
    # rpasblast calls here
    
    
    
}
