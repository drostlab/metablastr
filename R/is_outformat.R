is_outformat <- function(out.format) {
    if (!is.element(out.format, list_outformats()))
        stop("'",out.format, "' is not an available BLAST output format.", call. = FALSE)
    
    return(TRUE)
}
