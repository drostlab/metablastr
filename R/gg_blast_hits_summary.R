#' @title Plot \code{\link{blast_nucleotide_to_genomes}} Result
#' @description This function generates a joyplot of the \code{\link{blast_nucleotide_to_genomes}} output.  
#' @param blast_tbl a BLAST table generated with \code{\link{blast_nucleotide_to_genomes}}.
#' @param scope_cutoff The scope is defined as \code{1 - (abs(q_len - alig_length) / q_len))}. The \code{scope_cutoff}
#' defines the minimum scope that is required to retain a BLAST hit. Default is \code{scope_cutoff = 0.1} (meaning that each BLAST hit must have at least 0.1 scope). 
#' @param alpha a value passed to \code{\link[ggplot2]{aes}} specifying the transparency of distributions.
#' @param scale the ridges scale passed to \code{\link[ggridges]{geom_density_ridges}}.
#' @param query_name name of the query in the \code{blast_tbl}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title figure title which will be concatenated with the number of BLAST hits available in the inpit \code{blast_tbl}.
#' @param xticks number of ticks on the x-axis. Default is \code{xticks = 5}.
#' @param levels a character vector specifying the exact order of species names (levels) 
#' that is used to \code{\link{factor}} and sort species in the \code{\link[ggridges]{geom_density_ridges}} plot.
#' @author Hajk-Georg Drost
#' @export
gg_blast_hits_summary <- function(blast_tbl,
                                  scope_cutoff = 0,
                                  alpha = 0.7,
                                  scale = 4,
                                  query_name = "query",
                                  xlab = "Query",
                                  ylab = "Density over Number of BLAST Hits",
                                  title = "Total number of BLAST hits: ",
                                  xticks = 5,
                                  levels = NULL) {
  
  if (scope_cutoff > 1)
    stop("Please specify a scope_cutoff between [0,1].", call. = FALSE)
  
  p1 <-
    metablastr::gg_blast_hits(
      blast_tbl,
      type = "scope",
      scope_cutoff = scope_cutoff,
      levels = names(table(blast_tbl$species)),
      xlab = paste0("Length homology to ", query_name, " in %"),
      trim = FALSE
    )
  
  p2 <-
    metablastr::gg_blast_hits(
      blast_tbl,
      type = "alig_length",
      scope_cutoff = scope_cutoff,
      levels = names(table(blast_tbl$species)),
      xlab = paste0("Alignment length with ", query_name, " in number of amino acids"),
      title =  NULL,
      trim = FALSE
    )
  
  species <- scope <- NULL
  p3 <-
    ggplot2::ggplot(dplyr::filter(blast_tbl, scope >= scope_cutoff), ggplot2::aes(species)) + ggplot2::geom_bar() + ggplot2::coord_flip() + ggplot2::ylab("Number of homologous BLAST hits") +
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
    ) + ggplot2::xlab("")
  
  
  p_final <-
    gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
  
  return(p_final)
}
