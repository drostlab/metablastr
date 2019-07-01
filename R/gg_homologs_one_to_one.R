#' @title Plot number of BLAST hits per species as barplot from a \code{tibble} generated with \code{detect_homologs_genome_to_genome} or \code{detect_homologs_proteome_to_proteome}
#' @description This function generates a boxplot of different types of the \code{\link{detect_homologs_proteome_to_proteome}} output.
#' @param blast_tbl a BLAST table generated with \code{\link{detect_homologs_proteome_to_proteome}}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title plot title.
#' @param text_size size of label text.
#' @param y_ticks number of ticks on the y-axis.
#' @author Hajk-Georg Drost
#' @export

gg_homologs_one_to_one <- function(blast_tbl,
                                   xlab = "Species",
                                   ylab = "Number of one-to-one homologs (including splice variants)",
                                   title = "One-to-one homologs",
                                   text_size = 18,
                                   y_ticks = 8) {
  
  species <- query_id <- n_hits <- NULL
  
  blast_tbl_summary <- dplyr::summarize(dplyr::group_by(blast_tbl, species, query_id), n_hits = dplyr::n())
  blast_tbl_summary <- dplyr::filter(blast_tbl_summary, n_hits == 1)
  p <- ggplot2::ggplot(
    blast_tbl_summary,
    ggplot2::aes(x = species, y = n_hits, fill = species)
  ) +
    ggplot2::geom_bar(stat = "identity", show.legend = FALSE) +
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
    )) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks)) +
    ggsci::scale_fill_lancet() + ggplot2::ggtitle(title)
  
  return(p)
}
