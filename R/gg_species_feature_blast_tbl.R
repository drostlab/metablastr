#' @title Plot pairwise feature distributions of orthologous genes
#' @description This function enables users to plot pairwise feature distributions of orthologous genes.
#' @param blast_tbl a \code{blast_tbl}.
#' @param type type of sustitution rate quantification that shall be visualized on the y-axis. Options are:
#' \itemize{
#' \item \code{type = "q_len"}
#' \item \code{type = "qcov"}
#' \item \code{type = "qcovhsp"}
#' \item \code{type = "perc_identity"}
#' \item \code{type = "alig_length"}
#' \item \code{type = "mismatches"}
#' \item \code{type = "gap_openings"}
#' \item \code{type = "n_gaps"}
#' \item \code{type = "evalue"}
#' \item \code{type = "bit_score"}
#' \item \code{type = "score_raw"}
#' }
#' @param order a character vector specifying the order of the species using their scientific names (which have to match the species names stored in the \code{subject_species}
#' column of the input \code{blast_tbl}). 
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param title title label.
#' @author Hajk-Georg Drost
#' @export

gg_species_feature_blast_tbl <- function(blast_tbl, type, order, xlab = type, ylab = "Density", title = "") {
  
  blast_tbl_types <- c(
    "q_len",
    "qcov",
    "qcovhsp",
    "perc_identity",
    "alig_length",
    "mismatches",
    "gap_openings",
    "n_gaps",
    "s_len",
    "evalue",
    "bit_score",
    "score_raw"
  )
  
  if (!is.element(
    type,
    blast_tbl_types
  ))
  stop(
    "Please specify a type that is supported by this function: ", paste0(blast_tbl_types, collapse = ", "),
    call. = FALSE
  )
  
  subject_species <- q_len <- qcov <- qcovhsp <- perc_identity <- alig_length <- ..density.. <- NULL
  mismatches <- gap_openings <- n_gaps <- s_len <- evalue <- bit_score <- score_raw <- NULL
    
  if (type == "q_len") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = q_len,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "qcov") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = qcov,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "qcovhsp") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = qcovhsp,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "perc_identity") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = perc_identity,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "alig_length") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = alig_length,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "mismatches") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = mismatches,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "gap_openings") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = gap_openings,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "n_gaps") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = n_gaps,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "s_len") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = s_len,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "evalue") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = evalue,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "bit_score") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = bit_score,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  if (type == "score_raw") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = score_raw,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  
  p <- p +
    ggridges::geom_density_ridges(
      scale = 5,
      show.legend = FALSE,
      alpha = 0.4,
      size = 1.5,
      stat = "density", trim = TRUE
    ) + ggridges::theme_ridges() + ggsci::scale_colour_lancet() + ggsci::scale_fill_lancet() + ggplot2::labs(x = xlab,
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
    ) + ggplot2::ggtitle(title) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
  
  return(p)
}
