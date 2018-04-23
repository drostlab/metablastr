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
#' @param cores number of cores for parallel rpsblast searches.
#' @param blast.path path to rpsblast executables.
#' @author Hajk-Georg Drost, Anna Gogleva
#' @seealso \code{\link{blast_protein_to_protein}}, \code{\link{blast_nucleotide_to_protein}}
#' @export
#' @examples
#' small_fasta <- '/home/anna/anna/Labjournal/SecretSanta_external/test_fastas/small_10.fasta'
#' db_path <- '/home/anna/anna/Labjournal/db_for_protoro'
#' db_alias <- 'myCdd'
#' result <- blast_rpsblast(query = small_fasta, db = db_path, db.alias = db_alias)
#' 

query = small_fasta
db = db_path
db.alias = db_alias
output.path = NULL
db.import = FALSE
postgress.user = NULL
evalue = 1E-3
out.format = 8
cores = 1
blast.path = NULL

blast_rpsblast <- function(query,
                           db,
                           db.alias,
                           output.path = NULL,
                           db.import = FALSE,
                           postgress.user = NULL,
                           evalue = 1E-3,
                           out.format = 8,
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
    
    message("Starting 'rpsblast","' with  query: ", query, " against database: ", db.alias," using ", cores, " core(s) ...")
    
    # rpasblast calls here
    
    rpsblast_call <- paste('rpsblast -i',
                           query,
                           '-d', 
                           full.db.path
                           )
    # this does not work yet
    
    # output_rpsblast <-
    #     file.path(ifelse(is.null(output.path), ws_wrap(getwd()), ws_wrap(output.path)),
    #               paste0(unlist(stringr::str_split(
    #                   basename(query), "[.]"
    #               ))[1], ".blast_tbl"))
    
    # the call itself
    system(
        paste0(
            ifelse(is.null(blast.path), rpsblast_call, paste0("export PATH=$PATH:", rpsblast_call)),
               ' -e ',
               evalue,
               " -m 8 ",
               ' -a ',
               cores
               ))
        
}
