gg_shared_splice_variants <- function(shared_genes_tbl, 
                            xlab = "Splice Variants conserved (homologous) in at least N other species",
                            ylab = "Number of homologous splice variants",
                            title = "",
                            text_size = 18,
                            y_ticks = 8) {
  
  p <- ggplot2::ggplot(shared_genes_tbl, ggplot2::aes(x = factor(min_species, levels = 1:max(shared_genes_tbl$min_species)), y = n_splice_var)) +
    ggplot2::geom_bar(stat = "identity", fill = "#00A087B2") +
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