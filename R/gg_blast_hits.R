#' @title Plot \code{\link{blast_genomes}} Result
#' @description This function generates a joyplot of the \code{\link{blast_genomes}} output.  
#' @param blast_tbl a BLAST table generated with \code{\link{blast_genomes}}.
#' @param type the type of BLAST hit visualization. Options are:
#' \itemize{
#' \item \code{type = "scope"} 
#' \item \code{type = "alig_length"}
#' \item \code{type = "evalue"}
#' \item \code{type = "bit_score"}
#' \item \code{type = "qcovhsp"}
#' \item \code{type = "qcov"}
#' \item \code{type = "perc_identity"}
#' }
#' @param scope_cutoff The scope is defined as \code{1 - (abs(q_len - alig_length) / q_len))}. The \code{scope_cutoff}
#' defines the minimum scope that is required to retain a BLAST hit. Default is \code{scope_cutoff = 0.1} (meaning that each BLAST hit must have at least 0.1 scope). 
#' @param alpha a value passed to \code{\link[ggplot2]{aes}}specifying the transparency of distributions.
#' @param scale the ridges scale passed to \code{\link[ggridges]{geom_density_ridges}}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title figure title which will be concatenated with the number of BLAST hits available in the inpit \code{blast_tbl}.
#' @param xticks number of ticks on the x-axis. Default is \code{xticks = 5}.
#' @param levels a character vector specifying the exact order of species names (levels) 
#' that is used to \code{\link{factor}} and sort species in the \code{\link[ggridges]{geom_density_ridges}} plot.
#' @author Hajk-Georg Drost
#' @export
gg_blast_hits <-
  function(blast_tbl,
           type = "scope",
           scope_cutoff = 0.1,
           alpha = 0.7,
           scale = 4,
           xlab = "Query",
           ylab = "Density over Number of BLAST Hits",
           title = "Total number of BLAST hits: ",
           xticks = 5,
           levels = NULL
           ) {
    if (!is.element(type, c("scope", "alig_length", "evalue", "bit_score", "qcovhsp", "qcov", "perc_identity")))
      stop(
        "The argument '",
        type,
        "' is not a valid type option. Please consult the documentation for details."
      )
    
    q_len <- alig_length <- scope <- bit_score <- species <- evalue <- ..density.. <- NULL
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
            alpha = alpha,
            height = ..density..
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
            alpha = alpha,
            height = ..density..
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
            alpha = alpha,
            height = ..density..
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
            alpha = alpha,
            height = ..density..
          )
        )
    }
    
    if (type == "qcovhsp") {
      qcovhsp <- NULL
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = qcovhsp,
            y = species,
            fill = species,
            alpha = alpha,
            height = ..density..
          )
        )
    }
    
    if (type == "qcov") {
      qcovhsp <- NULL
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = qcov,
            y = species,
            fill = species,
            alpha = alpha,
            height = ..density..
          )
        )
    }
    
    if (type == "perc_identity") {
      qcovhsp <- NULL
      p <-
        ggplot2::ggplot(
          dplyr::filter(blast_tbl, scope >= scope_cutoff),
          ggplot2::aes(
            x = perc_identity,
            y = species,
            fill = species,
            alpha = alpha,
            height = ..density..
          )
        )
    }
    
    p <- p +
      ggridges::geom_density_ridges(scale = scale,
                                    show.legend = FALSE,
                                    alpha = alpha,
                                    size = 1.5,
                                    stat = "density", trim = TRUE) + ggridges::theme_ridges() +
      ggridges::theme_ridges(font_size = 13, grid = TRUE)  + 
      ggsci::scale_colour_lancet() + ggsci::scale_fill_lancet() +
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
      ) + ggplot2::ggtitle(paste0(title, nrow(
        dplyr::filter(blast_tbl, scope >= scope_cutoff)
      ))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 1,
        hjust = 1
      ))
      
    return(p)
  }
