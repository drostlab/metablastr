#' @title Plot number of pairwise orthologs for multiple species comparisons
#' @description Visualizes pairwise orthologs for multiple species comparisons as line graph.
#' @param ortho_summary a \code{data.frame} storing the \code{subject_species} in the first column and the number of orthologs named \code{n_orthologs} in the second column.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param title plot title.
#' @param vline plot a line characterizing the core set. Default is \code{vline = NULL} meaning that no vline is drawn.
#' @param ymax max value of y-axis. Default is \code{ymax = NULL} meaning that the largest value within the dataset will be used to determine y-max.
#' @author Hajk-Georg Drost
#' @export

gg_pairwise_orthologs_line <-
  function(ortho_summary,
           xlab = "Subject Species",
           ylab = "Number of reciprocal best hit orthologs",
           title = "",
           vline = NULL,
           ymax = NULL) {
    
    if (names(ortho_summary)[1] != "subject_species" || names(ortho_summary)[2] != "n_orthologs")
      stop("Please provide a ortho_summary data.frame or tibble with 2 columns and column names: 'subject_species' and 'n_orthologs'.", call. = FALSE)
    
    subject_species <- n_orthologs <- NULL
    
    p <- ggplot2::ggplot(ortho_summary,
                         ggplot2::aes(x = subject_species,
                                      y = n_orthologs,
                                      group = 1)) + ggplot2::geom_line(size = 2) + ggplot2::geom_point(size = 4)
      
      if (!is.null(vline)) {
        p <- p + ggplot2::geom_abline(
          intercept = vline,
          size = 2,
          col = "darkred",
          alpha = 0.4
        )
      }
      
    if (is.null(ymax))
      ymax <- round(max(ortho_summary$n_orthologs) + (max(ortho_summary$n_orthologs) * 0.15), 0)
    
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = n_orthologs),
        hjust = 0,
        vjust = -1.5,
        size = 6
      ) + ggplot2::scale_y_continuous(limits = c(0, ymax),
                                      breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = xlab, y = ylab, title = title) +
      ggplot2::theme(
        title            = ggplot2::element_text(size = 18, face = "bold"),
        legend.title     = ggplot2::element_text(size = 18, face = "bold"),
        legend.text      = ggplot2::element_text(size = 18, face = "bold"),
        axis.title       = ggplot2::element_text(size = 18, face = "bold"),
        axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
        axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
        panel.background = ggplot2::element_blank(),
        strip.text.x     = ggplot2::element_text(
          size           = 18,
          colour         = "black",
          face           = "bold"
        )
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 1,
        hjust = 1
      ))
    
    
    return(p)
  }
