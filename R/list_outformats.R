#' @title List available BLAST output formats
#' @description Available BLAST output formats are listed as they can
#' be specified as \code{out.format} argument in any \code{blast_*} function.
#' @author Hajk-Georg Drost
#' @examples 
#' # list available BLAST output formats
#' list_outformats()
#' @export
list_outformats <- function() {
    available_outformats <- c(
        "pair",
        "qa.ident",
        "qa.nonident",
        "fq.ident",
        "fq.nonident",
        "xml",
        "tab",
        "tab.comment",
        "ASN.1.text",
        "ASN.1.binary",
        "csv",
        "ASN.1",
        "json.seq.aln",
        "json.blast.multi",
        "xml2.blast.multi",
        "json.blast.single",
        "xml2.blast.single",
        "report"
    )
    
    return(available_outformats)
}
