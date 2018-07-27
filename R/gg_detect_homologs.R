#' @title 
#' @description 
#' @param 
#' @param 
#' @param 
gg_detect_homologs <- function(blast_tbl, type = "perc_identity", xlab = "Species", ylab = "", text_size = 18, y_ticks = 8) {
  
  if (!is.element(type, c("perc_identity", "alig_length", "scope")))
    stop("Please choose available types: type = 'perc_identity' ; type = 'alig_length' ; type = 'scope'.", call. = FALSE)
  
  if (type == "perc_identity")
  p <- ggplot2::ggplot(blast_tbl, ggplot2::aes(x = species, y = perc_identity, fill = species))
    
  if (type == "alig_length")
    p <- ggplot2::ggplot(blast_tbl, ggplot2::aes(x = species, y = alig_length, fill = species))
  
  if (type == "scope")
    p <- ggplot2::ggplot(blast_tbl, ggplot2::aes(x = species, y = scope, fill = species))
  
  p <- p +  
    ggplot2::geom_boxplot()  +
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
    ggsci::scale_fill_lancet()
  
  return(p)
}


