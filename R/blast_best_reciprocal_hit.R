#' @title Retrieve only the best reciprocal BLAST hit for each query
#' @description This function performs a BLAST search between query and subject sequences reciprocally and
#' returns (for each direction) only the best hit based on the following criteria.
#' 
#' A best blast hit is defined as:
#' \itemize{
#' \item the hit with the smallest e-value
#' \item if e-values are identical then the hit with the longest alignment length is chosen
#' }
#' 
#' A best reciprocal BLAST hit is defined as a hit that fulfills the following symmetry criteria:
#' 
#' \code{blast_best_hit(A,B) = blast_best_hit(B,A)}
#' 
#' @param query path to input file in fasta format.
#' @param subject path to subject file in fasta format or blast-able database.
#' @param search_type type of query and subject sequences that will be compared via BLAST search.
#' Options are:
#' \itemize{
#' \item \code{search_type = "nucleotide_to_nucleotide"}
#' \item \code{search_type = "nucleotide_to_protein"}
#' \item \code{search_type = "protein_to_nucleotide"}
#' \item \code{search_type = "protein_to_protein"}
#' }
#' @param strand Query DNA strand(s) to search against database/subject.
#' Options are:
#' \itemize{
#' \item \code{strand = "both"} (Default): query against both DNA strands.
#' \item \code{strand = "minus"} : query against minus DNA strand.
#' \item \code{strand = "plus"} : query against plus DNA strand.
#' }
#' @param output.path path to folder at which BLAST output table shall be stored. 
#' Default is \code{output.path = NULL} (hence \code{getwd()} is used).
#' @param is.subject.db logical specifying whether or not the \code{subject} file is a file in fasta format (\code{is.subject.db = FALSE}; default)
#' or a blast-able database that was formatted with \code{makeblastdb} (\code{is.subject.db = TRUE}).
#' @param task BLAST search task option (depending on the selected \code{search_type}). Options are:
#' \itemize{
#' \item \code{search_type = "nucleotide_to_nucleotide"}
#' \itemize{
#' \item \code{task = "blastn"} : Standard nucleotide-nucleotide comparisons (default) - Traditional BLASTN requiring an exact match of 11.
#' \item \code{task = "blastn-short"} : Optimized nucleotide-nucleotide comparisons for query sequences shorter than 50 nucleotides.
#' \item \code{task = "dc-megablast"} : Discontiguous megablast used to find somewhat distant sequences.
#' \item \code{task = "megablast"} : Traditional megablast used to find very similar (e.g., intraspecies or closely related species) sequences.
#' \item \code{task = "rmblastn"}
#' }
#' \item \code{search_type = "nucleotide_to_protein"}
#' \itemize{
#' \item \code{task = "blastx"} : Standard nucleotide-protein comparisons (default).
#' \item \code{task = "blastx-fast"} : Optimized nucleotide-protein comparisons.
#' }
#' \item \code{search_type = "protein_to_nucleotide"}
#' \itemize{
#' \item \code{task = "tblastn"} : Standard protein-nucleotide comparisons (default).
#' \item \code{task = "tblastn-fast"} : Optimized protein-nucleotide comparisons.
#' }
#' \item \code{search_type = "protein_to_protein"}
#' \itemize{
#' \item \code{task = "blastp"} : Standard protein-protein comparisons (default).
#' \item \code{task = "blast-fast"} : Improved BLAST searches using longer words for protein seeding.
#' \item \code{task = "blastp-short"} : Optimized protein-protein comparisons for query sequences shorter than 30 residues.
#' }
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
#' test_best_reciprocal_hit <- blast_best_reciprocal_hit(
#'                  query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
#'                  subject = system.file('seqs/sbj_nn_best_hit.fa', package = 'metablastr'),
#'                  search_type = "nucleotide_to_nucleotide",
#'                  output.path = tempdir(),
#'                  db.import  = FALSE)
#'                  
#'  # look at results
#'  test_best_reciprocal_hit
#' }
#' 
#' @seealso \code{\link{blast_nucleotide_to_nucleotide}}, \code{\link{blast_protein_to_protein}},
#' \code{\link{blast_nucleotide_to_protein}}, \code{\link{blast_protein_to_nucleotide}}, \code{\link{blast_best_hit}}
#' @export

blast_best_reciprocal_hit <-
  function(query,
           subject,
           search_type = "nucleotide_to_nucleotide",
           strand = "both",
           output.path = NULL,
           is.subject.db = FALSE,
           task = "blastn",
           db.import = FALSE,
           postgres.user = NULL,
           evalue = 0.001,
           out.format = "csv",
           cores = 1,
           max.target.seqs = 10000,
           db.soft.mask = FALSE,
           db.hard.mask = FALSE,
           blast.path = NULL) {
    
  
  
  if (search_type == "nucleotide_to_nucleotide") {
    
    orthoA <-
      blast_best_hit(
        query = query,
        subject = subject,
        search_type = "nucleotide_to_nucleotide",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
    orthoB <-
      blast_best_hit(
        query = subject,
        subject = query,
        search_type = "nucleotide_to_nucleotide",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
  }
  
  if (search_type == "nucleotide_to_protein") {
    
    orthoA <-
      blast_best_hit(
        query = query,
        subject = subject,
        search_type = "nucleotide_to_protein",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
    orthoB <-
      blast_best_hit(
        query = subject,
        subject = query,
        search_type = "protein_to_nucleotide",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
  }
  
  if (search_type == "protein_to_nucleotide") {
    
    orthoA <-
      blast_best_hit(
        query = query,
        subject = subject,
        search_type = "protein_to_nucleotide",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
    orthoB <-
      blast_best_hit(
        query = subject,
        subject = query,
        search_type = "nucleotide_to_protein",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
  }
  
  if (search_type == "protein_to_protein") {
    
    orthoA <-
      blast_best_hit(
        query = query,
        subject = subject,
        search_type = "protein_to_protein",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
    orthoB <-
      blast_best_hit(
        query = subject,
        subject = query,
        search_type = "protein_to_protein",
        strand = strand,
        output.path = output.path,
        is.subject.db = is.subject.db,
        task = task,
        db.import = db.import,
        postgres.user = postgres.user,
        evalue = evalue,
        out.format = out.format,
        cores = cores,
        max.target.seqs = max.target.seqs,
        db.soft.mask = db.soft.mask,
        db.hard.mask = db.hard.mask,
        blast.path = blast.path
      )
    
  }
  
  colnames(orthoB)[1:2] <- c("subject_id", "query_id")
  
  tryCatch({
    return(dplyr::semi_join(
      orthoA,
      orthoB,
      by = c("query_id", "subject_id")
    ))
    
    
  }, error = function(e) {
    stop(
      "The BLAST tables resulting from ",
      query,
      " and ",
      subject,
      " could not be joined properly to select only the reciprocal best hits."
    )
  })
  
  
  
  
}
