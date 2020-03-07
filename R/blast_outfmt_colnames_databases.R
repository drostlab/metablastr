blast_outfmt_colnames_databases <- function() {
  
  params <-
    c(
      "query_id",
      "subject_id",
      "perc_identity",
      "num_ident_matches",
      "alig_length",
      "mismatches",
      "gap_openings",
      "n_gaps",
      "pos_match",
      "ppos",
      "q_start",
      "q_end",
      "q_len",
      "qcov",
      "qcovhsp",
      "s_start",
      "s_end",
      "s_len",
      "evalue",
      "bit_score",
      "score_raw",
      "sscinames",
      "staxids",
      "scomnames",
      "sblastnames",
      "sskingdoms",
      "sstrand"
    )
  
  return(params)
}
