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
#' @param out.format a character string specifying the format of the file in which the rpsblast results shall be stored.
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
#'  }
#' @param cores number of cores for parallel rpsblast searches.
#' @param blast.path path to rpsblast executables.
#' @author Hajk-Georg Drost, Anna Gogleva
#' @seealso \code{\link{blast_protein_to_protein}}, \code{\link{blast_nucleotide_to_protein}}
#' @export
#' @examples
#' db_path <- '/home/anna/anna/Labjournal/db_for_protoro'
#' db_alias <- 'myCdd'
#' rps_test <- blast_rpsblast(query = system.file('seqs/sbj_aa.fa', package = 'metablastr'), db = db_path, db.alias = db_alias, output.path = tempdir())

query = system.file('seqs/sbj_aa.fa', package = 'metablastr')

blast_rpsblast <- function(query,
                           db,
                           db.alias,
                           output.path = NULL,
                           db.import = FALSE,
                           postgres.user = NULL,
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
    db.postfix <- c('aux', 'freq', 
                    'loo', 'phr',
                    'pin', 'psi',
                    'psq', 'rps')

    full.db.path <- paste(db, db.alias, sep = '/')
    db.files <- sapply(db.postfix,
                       function(x) paste(full.db.path, x, sep = '.'))
    db.status <- sapply(db.files, file.exists)
    
    if (!(all(db.status))) {
        missings <- paste(names(db.status[db.status == FALSE]), collapse = ' ')
        stop('Unfortunately, one ore more required database files are missing: ', missings, call. = FALSE)
    }
    
    message("Starting 'rpsblast+","' with  query: ", query, " against database: ", db.alias," using ", cores, " core(s) ...")
    
    # rpasblast call starts here
    
    rpsblast_call <- paste('rpsblast+ -query',
                           query,
                           '-db', 
                           full.db.path
                           )
  
    # output file with ws handling
    
    output_rpsblast <-
        file.path(ifelse(is.null(output.path), ws_wrap(getwd()), ws_wrap(output.path)),
                  paste0(unlist(stringr::str_split(
                      basename(query), "[.]"
                  ))[1], ".rpsblast_tbl"))
    
    # don't need ws handling anymore when reading file
    output_read_rpsblast <-
        file.path(ifelse(is.null(output.path), getwd(), output.path),
                  paste0(unlist(stringr::str_split(
                      basename(query), "[.]"
                  ))[1], ".rpsblast_tbl"))
    
    # prep database if not ready yet (to do)
    
    
    # the rpsblast call itself
    system(
        paste0(
            ifelse(is.null(blast.path), rpsblast_call, paste0("export PATH=$PATH:", rpsblast_call)),
               ' -evalue ',
               evalue,
               paste0(' -outfmt ', outformat2num(out.format = out.format)),
               ' -num_threads ',
               cores,
               ' -out ',
               output_rpsblast
               ))
    
    # read rpasblast outputs and convert into tibble (to do)
}