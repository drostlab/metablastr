gg_shared_genes <- function(shared_genes_tbl, 
                            xlab = "Genes conserved (homologous) in at least N other species",
                            ylab = "Number of homologous genes",
                            title = "",
                            text_size = 18,
                            y_ticks = 8) {
  
  min_species <- n_genes <- NULL
  p <- ggplot2::ggplot(shared_genes_tbl, ggplot2::aes(x = factor(min_species, levels = 1:max(shared_genes_tbl$min_species)), y = n_genes)) +
    ggplot2::geom_bar(stat = "identity", fill = "#2D2D2D") +
    ggplot2::facet_grid(. ~ qcovhsp) +
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
    ) + ggplot2::ggtitle(title) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks)) 
  
  return(p)
}
