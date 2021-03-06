% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast_protein_to_proteomes.R
\name{blast_protein_to_proteomes}
\alias{blast_protein_to_proteomes}
\title{Perfrom BLAST Searches Against a Set of Proteomes}
\usage{
blast_protein_to_proteomes(
  query,
  subject_proteomes,
  blast_type = "blastp",
  blast_output_path = "blast_output",
  min_alig_length = 30,
  evalue = 1e-05,
  max.target.seqs = 5000,
  update = FALSE,
  ...
)
}
\arguments{
\item{query}{path to input file in fasta format.}

\item{subject_proteomes}{a vector containing file paths to the reference proteomes that shall be queried (e.g. file paths returned by \code{\link[biomartr]{meta.retrieval}}).}

\item{blast_type}{specification of the BLAST type shall be used to perform BLAST searches between query and reference.
Available options are:
\itemize{
\item \code{task = "blastp"} : Standard protein-protein comparisons (default); (see \code{\link{blast_protein_to_protein}} for details).
\item \code{task = "blast-fast"} : Improved BLAST searches using longer words for protein seeding.
\item \code{task = "blastp-short"} : Optimized protein-protein comparisons for query sequences shorter than 30 residues.
}}

\item{blast_output_path}{a path to a folder that will be created to store BLAST output tables for each individual query-proteome search.}

\item{min_alig_length}{minimum alignment length that shall be retained in the result dataset. All hit alignments with smaller
hit alignment length will be removed automatically.}

\item{evalue}{minimum expectation value (E) threshold for retaining hits (default: evalue = 0.00001).}

\item{max.target.seqs}{maximum number of aligned sequences that shall be retained. Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .}

\item{update}{a logical value indicating whether or not pre-computed BLAST tables should be removed and re-computed (\code{update = TRUE})
or imported from existing file (\code{update = FALSE}) (Default).}

\item{\dots}{additional arguments passed to \code{\link{blast_protein_to_protein}}.}
}
\description{
This function takes a fasta file containing query sequences
as input and performs BLAST searches of these query sequences against
a set of reference proteomes to retrieve corresponding hits in diverse proteomes
}
\details{
The \code{blast_protein_to_proteomes} function enables users to BLAST specific query sequences against a set of reference proteomes
and retrieve the corresponding BLAST output.
}
\author{
Hajk-Georg Drost
}
