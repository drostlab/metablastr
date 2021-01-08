# NEWS

## metablastr v0.3.1

- fixing bug in `extract_promotor_seqs_from_genome()` which causes genes on the minus strand that have `end - start + 1 > 0` coordinates in their respective annotation file to be removed from the analysis. The new convention is now that genes will be filtered according to the `strand` column
of the annotation file (hence: `+` or `-`) and it will be assumed that both, plus and minus strands use `start` and `end` coordinates according to their increasing direction, e.g. `start = 1` and `end = 5000`. (Many thanks for detecting this: Alexander Gabel) 