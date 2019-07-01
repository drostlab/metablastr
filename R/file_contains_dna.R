#' @title Helper function to check whether input file contains dna as required
#' @description Check whether input file contains dna as required.
#' @param file file path to \code{fasta} file.
#' @param input_type a character string specifying whether input file a \code{"query"} or \code{"subject"} file.
#' @author Hajk-Georg Drost
#' @return When a file contains the correct sequence type then this function
#' will only execute without any output. In case a file does not contains the
#' correct sequence type then an error is thrown.
#' @examples
#' # test an dna file
#' file_contains_dna(file = system.file('seqs/qry_nn.fa', package = 'metablastr'),
#'                   input_type = "query")
#' @export

file_contains_dna <- function(file, input_type) {
  if (!is.element(input_type, c("query", "subject")))
    stop("Please provide a valid input_type.", call. = FALSE)
  
  alphabet <-
    paste0(unique(stringr::str_to_upper(Biostrings::uniqueLetters(
      Biostrings::readBStringSet(filepath = file, nrec = 1)
    ))), collapse = "")
  
    if (!any(is.element(alphabet, c("ACGNT", "ACGT")))) {
      message("Your ",
              input_type,
              " file contains the folowing letters: '",
              alphabet,
              "'.")
      message("However, these letters shoould be present: ", "'ACGT' or 'ACGNT'", ".")
      stop(
        "Please provide the correct sequence type for the ",
        input_type,
        " file. The correct sequence type is: DNA.",
        call. = FALSE
      )
    }
}
