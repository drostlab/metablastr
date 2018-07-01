gg_blast_hits <-
  function(blast_tbl,
           type = "scope",
           scope_cutoff = 0.1,
           alpha = 0.7,
           scale = 4,
           xlab = "Query",
           ylab = "Density over Number of BLAST Hits",
           xticks = 5,
           levels = NULL
           ) {
    if (!is.element(type, c("scope", "alig_length", "evalue", "bit_score")))
      stop(
        "The argument '",
        type,
        "' is not a valid type option. Please consult the documentation for details."
      )
    
    q_len <- alig_length <- scope <- NULL
    blast_tbl <-
      dplyr::mutate(blast_tbl, scope = 1 - (abs(q_len - alig_length) / q_len))
    
    if (!is.null(levels))
      blast_tbl <- dplyr::mutate(blast_tbl, species = factor(species, levels = levels))
    
    if (type == "scope") {
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = scope * 100,
            y = species,
            fill = species,
            alpha = alpha
          )
        )
    }
    
    if (type == "alig_length") {
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = alig_length,
            y = species,
            fill = species,
            alpha = alpha
          )
        )
    }
    
    if (type == "evalue") {
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = evalue,
            y = species,
            fill = species,
            alpha = alpha
          )
        )
    }
    
    if (type == "bit_score") {
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = bit_score,
            y = species,
            fill = species,
            alpha = alpha
          )
        )
    }
    
    p <- p +
      ggridges::geom_density_ridges(scale = scale, show.legend = FALSE) + ggridges::theme_ridges() +
      ggridges::theme_ridges(font_size = 13, grid = TRUE)  +
      ggplot2::labs(x = xlab,
                    y = ylab) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 18)) +
      ggplot2::theme(
        axis.title  = ggplot2::element_text(size = 18, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 18, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 18, face = "bold"),
        panel.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          size = 18,
          colour = "black",
          face = "bold"
        )
      ) + ggplot2::ggtitle(paste0("Total number of BLAST hits: ", nrow(
        dplyr::filter(blast_tbl, scope >= scope_cutoff)
      ))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 1,
        hjust = 1
      ))
      
    
    return(p)
  }