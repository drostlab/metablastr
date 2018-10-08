#' @title Boxplot of a \code{blast_tbl} generated with \code{detect_homologs_genome_to_genome} or \code{detect_homologs_proteome_to_proteome}
#' @description This function generates a boxplot of different types of the \code{\link{detect_homologs_proteome_to_proteome}} output.
#' @param blast_tbl a BLAST table generated with \code{\link{detect_homologs_proteome_to_proteome}}.
#' @param x_type variable thah shall be visualized on the x-axis. E.g.:
#' \itemize{
#' \item \code{type = "species"}: the species in which blast hits were found (Default)
#' }
#' @param y_type variable thah shall be visualized on the y-axis. E.g.:
#' \itemize{
#' \item \code{type = "qcovhsp"}: query coverage of the respective blast hit (Default)
#' \item \code{type = "perc_identity"}: visualize the alignment identity in percent for each BLAST hit.
#' \item \code{type = "alig_length"}: visualize the alignment length for each BLAST hit.
#' \item \code{type = "scope"}: visualize the length homology to the initial query for each BLAST hit.
#' }
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title plot title.
#' @param text_size size of label text.
#' @param y_ticks number of ticks on the y-axis.
#' @author Hajk-Georg Drost
#' @export

gg_hits_boxplot <-
  function(blast_tbl,
           x_type = "species",
           y_type = "qcovhsp",
           xlab = x_type,
           ylab = y_type,
           title = "",
           text_size = 18,
           y_ticks = 8) {
  
    
    if (!is.element(x_type, names(blast_tbl)))
      stop("Please choose a 'x_type' that is available in the 'blast_tbl': ",names(blast_tbl), call. = FALSE)
    
    if (!is.element(y_type, names(blast_tbl)))
      stop("Please choose a 'y_type' that is available in the 'blast_tbl': ",names(blast_tbl), call. = FALSE)
    
    
  p <- ggplot2::ggplot(blast_tbl, ggplot2::aes_string(x = x_type, y = y_type, colour = x_type))
    
  p <- p +  
    ggplot2::geom_boxplot(size = 1.5, alpha = 0.7)  +
    ggplot2::geom_point(ggplot2::aes_string(colour = x_type), alpha = 0.3) +
    ggplot2::geom_jitter(ggplot2::aes_string(colour = x_type), position = ggplot2::position_jitter(0.25), alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = xlab,
                  y = ylab) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = text_size)) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = text_size, face = "bold"),
      axis.text.y = ggplot2::element_text(size = text_size, face = "bold"),
      axis.text.x = ggplot2::element_text(size = text_size, face = "bold"),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = text_size,
        colour = "black",
        face = "bold"
      )
    ) + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 1,
      hjust = 1
    )) + ggplot2::ggtitle(paste0(title, " (total: ",nrow(blast_tbl)," hits)")) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks)) +
    ggsci::scale_colour_lancet()
  
  return(p)
}


