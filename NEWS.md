# NEWS

## metablastr v0.3.2

- updating to GPL-2 license 

## metablastr v0.3.1

- fixing bug in `extract_promotor_seqs_from_genome()` which causes genes on the minus strand that have `end - start + 1 > 0` coordinates in their respective annotation file to be removed from the analysis. The new convention is now that genes will be filtered according to the `strand` column
of the annotation file (hence: `+` or `-`) and it will be assumed that both, plus and minus strands use `start` and `end` coordinates according to their increasing direction, e.g. `start = 1` and `end = 5000`. (Many thanks for detecting this: Alexander Gabel) 
- passing through `blast.path` to internal check function `is_blast_installed()`
- when performing `blast_best_reciprocal_hit()`, now users can specify tasks in both directions when selecting `search_type = protein_to_nucleotide` (`task = c("tblastn", "blastx")`) or `search_type = nucleotide_to_protein` (`task = c("blastx", "tblastn")`) #10. Many thanks to @wcfung14.
```r
blast_test_reciprocal <- blast_best_reciprocal_hit(
    query   = system.file('seqs/qry_aa.fa', package = 'metablastr'), # protein sequence
    subject = system.file('seqs/sbj_nn_best_hit.fa', package = 'metablastr'), # nucleotide sequence
    search_type = "protein_to_nucleotide",
    task = c("tblastn", "blastx"),
    evalue = 0.000001,
    output.path = tempdir(),
    db.import  = FALSE)
```