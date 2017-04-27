#' @title Import BLAST output file
#' @description When performing BLAST searches with the \code{blast_*()} functions,
#' the corresponding BLAST output file can be imported into the current R session using this function.
#' All output formats given by BLAST are supported (see e.g. description for details).
#' @param file path to BLAST output file.
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
#'  @param spark_version specify Spark version that shall be installed locally. Default is \code{spark_version = "2.0.1"}.
#'  @author Hajk-Georg Drost
#'  @seealso \code{\link{blast_protein}}
#'  @importFrom RPostgreSQL
#'  @export
read_blast <- function(file, out.format, postgres.user = NULL) {

   if (out.format == "tab") {
       
       # # install local version of Spark if not available yet
       # if (nrow(sparklyr::spark_installed_versions()) == 0) {
       #     sparklyr::spark_install(version = spark_version)
       # }
       # 
       # # open local connection to spark
       # sparkconnect <- sparklyr::spark_connect(master = "local")
       # 
       # import_blast_tbl <- sparklyr::copy_to(sparkconnect, iris, 
       #                                       "spark_blast_tbl",
       #                                        overwrite = TRUE) 
       
       
       if (is.null(postgres.user))
           stop("Please specify a 'postgres.user' to import BLAST output into PostgresSQL database.")
       
       require(RPostgreSQL)
     
       postgres_filename <- paste0(unlist(stringr::str_split(basename(file),"[.]"))[1],"_postgres")
       
       connect_db <-
           DBI::dbConnect(
               DBI::dbDriver("PostgreSQL"),
               user = postgres.user,
               password = "",
               host = "localhost",
               port = 5432, 
               dbname = postgres.user)
       
       DBI::dbWriteTable(
           connect_db,
           name      = postgres_filename,
           value     = file,
           row.names = FALSE,
           header    = FALSE,
           sep       = "\t",
           overwrite = TRUE
       )
       
       blast_sql_db <-
           dplyr::src_postgres(
               dbname = postgres.user,
               host = "localhost",
               port = 5432,
               user = postgres.user,
               password = ""
           )
       
       blast_postgres <-
           dplyr::tbl(blast_sql_db, postgres_filename)
       
       on.exit({
           #sparklyr::spark_disconnect(sparkconnect)
           DBI::dbDisconnect(connect_db)
           
       })
       
       return(blast_postgres)
   }     
    
   
    if (out.format == "csv") {
        blast_csv <- readr::read_delim(file, delim = ",", 
                                       col_names = FALSE,
                                       col_types = readr::cols(
                                           "X1" = readr::col_character(),
                                           "X2" = readr::col_character(),
                                           "X3" = readr::col_double(),
                                           "X4" = readr::col_integer(),
                                           "X5" = readr::col_integer(),
                                           "X6" = readr::col_integer(),
                                           "X7" = readr::col_integer(),
                                           "X8" = readr::col_integer(),
                                           "X9" = readr::col_integer(),
                                           "X10" = readr::col_double(),
                                           "X11" = readr::col_integer(),
                                           "X12" = readr::col_integer(),
                                           "X13" = readr::col_integer(),
                                           "X14" = readr::col_integer(),
                                           "X15" = readr::col_integer(),
                                           "X16" = readr::col_integer(),
                                           "X17" = readr::col_double(),
                                           "X18" = readr::col_double(),
                                           "X19" = readr::col_double()
                                       ))
        colnames(blast_csv) <- blast_outfmt_colnames()
        return(blast_csv)
    }    
    
}
