#' @title Perform Nucleotide to Protein BLAST Searches (BLASTX)
#' @description Run nucleotide to protein BLAST of reference sequences
#' against a blast-able database or fasta file. Internally BLAST translates the 
#' nucleotide sequence into a protein sequence and then searches for hits. 
#' @param query path to input file in fasta format.
#' @param subject path to subject file in fasta format or blast-able database.
#' @param strand Query DNA strand(s) to search against database/subject.
#' Options are:
#' \itemize{
#' \item \code{strand = "both"} : query against both DNA strands.
#' \item \code{strand = "minus"} : query against minus DNA strand.
#' \item \code{strand = "plus"} : query against plus DNA strand.
#' }
#' @param output.path path to folder at which BLAST output table shall be stored. 
#' Default is \code{output.path = NULL} (hence \code{getwd()} is used).
#' @param is.subject.db logical specifying whether or not the \code{subject} file is a file in fasta format (\code{is.subject.db = FALSE}; default)
#' or a blast-able database that was formatted with \code{makeblastdb} (\code{is.subject.db = TRUE}).
#' @param task nucleotide search task option. Options are:
#' \itemize{
#' \item \code{task = "blastx"} : Standard nucleotide-protein comparisons (default).
#' \item \code{task = "blastx-fast"} : Optimized nucleotide-protein comparisons.
#' }
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
#' @param cores number of cores for parallel BLAST searches.
#' @param max.target.seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param db.soft.mask shall low complexity regions be soft masked? Default is \code{db.soft.mask = FALSE}.
#' @param db.hard.mask shall low complexity regions be hard masked? Default is \code{db.hard.mask = FALSE}.
#' @param blast.path path to BLAST executables.
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' blast_test <- blast_nucleotide_to_protein(
#'                  query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
#'                  subject = system.file('seqs/sbj_aa.fa', package = 'metablastr'),
#'                  output.path = tempdir(),
#'                  db.import  = FALSE)
#'                  
#'  # look at results
#'  blast_test
#' }
#' 
#' @seealso \code{\link{blast_nucleotide_to_nucleotide}}, \code{\link{blast_protein_to_protein}}, 
#' \code{\link{blast_protein_to_nucleotide}}, \code{\link{blast_best_hit}}
#' @export

blast_nucleotide_to_protein <- function(query, 
                             subject,
                             strand = "both",
                             output.path = NULL,
                             is.subject.db = FALSE,
                             task = "blastx",
                             db.import = FALSE,
                             postgres.user = NULL,
                             evalue   = 1E-3,
                             out.format = "csv", 
                             cores = 1,
                             max.target.seqs = 10000,
                             db.soft.mask = FALSE,
                             db.hard.mask = FALSE,
                             blast.path = NULL) {
    
    if (!is_blast_installed())
        stop("Please install a valid version of BLAST.", call. = FALSE)
    
    if (db.import) {
        if (!is.element(out.format, c("xml", "tab", "csv")))
            stop("Only output formats: 'xml', 'tab', or 'csv' can be imported.", call. = FALSE)
    }
    
    if (!is.element(strand, c("both", "plus", "minus")))
        stop("Please specify a strand option that is supported by BLAST: strand = 'both', strand = 'plus', strand = 'minus'.")
    
    # make sure that input files contain the correct sequence type    
    file_contains_dna(query, "query")  
    file_contains_aa(subject, "subject")
  
    # determine the number of cores on a multicore machine
    multi.cores <- parallel::detectCores()
    
    # in case one tries to use more cores than are available
    if (cores > multi.cores)
        stop("You chose more cores than are available on your machine.", call. = FALSE)
    
    # test if query file exists
    if (!file.exists(query))
        stop("Unfortunately, no query file has been found at ", query, call. = FALSE)
    
    # test if subject file or database exists
    if (!file.exists(subject))
        stop("Unfortunately, no subject file has been found at ", subject, call. = FALSE)
    
    if (!is.element(task, c("blastx", "blastx-fast")))
        stop("Please choose a nucleotide-protein comparison task that is supported by BLAST: task = 'blastx' or task = 'blastx-fast'.", call. = FALSE)
    
    message("Starting 'blastx -task ", task,"' with  query: ", query, " and subject: ",subject," using ", cores, " core(s) ...")
    
    blast_call <-
        paste0("blastx -query ", ws_wrap(query), " -db ", ws_wrap(subject))
    
    output_blast <-
      file.path(ifelse(is.null(output.path), ws_wrap(getwd()), ws_wrap(output.path)),
                paste0(unlist(stringr::str_split(
                  basename(query), "[.]"
                ))[1],"_",unlist(stringr::str_split(
                  basename(subject), "[.]"
                ))[1],"_",task,"_eval_",evalue, ".blast_tbl"))
    
    # output_blast without ws_wrap()
    output_read_blast <-
      file.path(ifelse(is.null(output.path), getwd(), output.path),
                paste0(unlist(stringr::str_split(
                  basename(query), "[.]"
                ))[1],"_",unlist(stringr::str_split(
                  basename(subject), "[.]"
                ))[1],"_",task,"_eval_",evalue, ".blast_tbl"))
    
    
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
                    " -input_type fasta -dbtype nucl -hash_index"
                )
            )
        }
    } 
    
    system(
        paste0(
            ifelse(is.null(blast.path), blast_call, paste0("export PATH=$PATH:", blast_call)),
            " -evalue ",
            evalue,
            " -strand ", 
            strand,
            " -max_target_seqs ",
            max.target.seqs,
            " -out ",
            output_blast ,
            " -num_threads ",
            cores,
            ifelse(db.soft.mask, " -db_soft_mask", ""),
            ifelse(db.hard.mask, " -db_hard_mask", ""),
            paste0( " -task ", task),
            paste0(' -outfmt "', outformat2num(out.format = out.format), ' qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp sstart send slen evalue bitscore score"')
        )
    )
    
    if (!is.null(db.import)) {
      if (db.import) {
        blast_tbl <- read_blast(file = output_read_blast, 
                                out.format = "postgres",
                                postgres.user = postgres.user)
        message("\n")
        message("BLAST search finished! A PostgreSQL  database connection to the BLAST output file has been generated. The BLAST output file can be found at: ", output_blast)
        return(blast_tbl)
      } else {
        blast_tbl <- read_blast(file = output_read_blast, 
                                out.format = out.format,
                                postgres.user = NULL)
        
        message("\n")
        message("BLAST search finished! The BLAST output file was imported into the running R session. The BLAST output file has been stored at: ", output_blast)
        return(blast_tbl)
      }
    } else {
      message("\n")
      message("BLAST search finished! BLAST output file has been stored at: ", output_blast)
    }
    
}

