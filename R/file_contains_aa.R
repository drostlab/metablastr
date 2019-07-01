#' @title Helper function to check whether input file contains amino acids as required
#' @description Check whether input file contains amino acids as required.
#' @param file file path to \code{fasta} file.
#' @param input_type a character string specifying whether input file a \code{"query"} or \code{"subject"} file.
#' @author Hajk-Georg Drost
#' @return When a file contains the correct sequence type then this function
#' will only execute without any output. In case a file does not contains the 
#' correct sequence type then an error is thrown.
#' @examples 
#' # test an aa file
#' file_contains_aa(file = system.file('seqs/qry_aa.fa', package = 'metablastr'), input_type = "query")
#' @export
file_contains_aa <- function(file, input_type) {
  if (!is.element(input_type, c("query", "subject")))
    stop("Please provide a valid input_type.", call. = FALSE)
  
  alphabet <-
    unique(stringr::str_to_upper(Biostrings::uniqueLetters(
      Biostrings::readBStringSet(filepath = file, nrec = 1)
    )))
  
  if (any((!is.element(stringr::str_to_upper(alphabet), Biostrings::AA_ALPHABET)))) {
    message("Your ",
            input_type,
            " file contains the folowing letters: '",
            alphabet,
            "'.")
    message("However, these letters shoould be present: ", Biostrings::AA_ALPHABET, ".")
    stop("Please provide the correct sequence type for the ", input_type, " file. The correct sequence type is: AMINO ACIDS.", call. = FALSE)
  }
    
  
  if (any(!is.element(c("Q", "E"), stringr::str_to_upper(alphabet)))) {
    message("Your ",
            input_type,
            " file contains the folowing letters: '",
            alphabet,
            "'.")
    message("However, these letters shoould be present: ", Biostrings::AA_ALPHABET, ", Q, E.")
    stop("Please provide the correct sequence type for the ", input_type, " file. The correct sequence type is: AMINO ACIDS.", call. = FALSE)
  }
}







