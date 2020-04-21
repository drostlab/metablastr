# metablastr

## An easy-to-use framework to perform massive sequence searches with R

### Motivation 

The exponentially growing number of available sequences in biological databases
revolutionizes the way modern life science research is conducted. Approximately
hundred thousand genomic sequences spanning diverse species from the tree of life
are currently publically available and free to access. It is now possible to
access and retrieve this data automatically using the R package [biomartr](https://github.com/ropensci/biomartr)
and the next step is to harness this wealth of sequence diversity to explore
and detect novel patterns of evolvability, variation, and disease emergence.

The R package [biomartr](https://github.com/ropensci/biomartr)
__solves the problem of retrieving this vast amount of biological sequence data__ in a standardized and computationally reproducible way and
the `metablastr` package aims to __solve the problem of performing massive scale
sequence searches__ in a standardized and computationally reproducible way. 

Both packages, `biomartr` and `metablastr` are designed to complement
each other seamlessly to provide users with a toolset to automatically
retrieve thousands of biological sequences (thousands of genomes, proteomes, annotations, etc)
and to use these sequences to perform massive sequence searches to
extract novel patterns of similarity and divergence between large sets
of species.

The most prominent tool to perform sequence searches at scale is the Basic Local Alignment Search Tool (BLAST)
which is designed to find regions of sequence similarity between query and a subject sequences or sequence databases.
However, in addition to BLAST recent efforts in bioinformatics research led to the
development of faster and/or more accurate sequence search tools which enable
sequence searches and sequence comparisons between thousands of genomes.

Examples include:

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [MMSeqs2](https://github.com/soedinglab/MMseqs2)
- [HH-suite3](https://github.com/soedinglab/hh-suite)
- [LAST+](https://github.com/hallamlab/LAST-Plus)

### Short package description  

The `metablastr` package harnesses the power of these search tools by providing interface functions between R and the standalone (command line tool) versions
of these programs. In addition to providing interface functions, `metablastr` provides a scalable database backend infrastructure
and analytics tools to store and handle the extensive search output generated
by these tools when handling thousands of genomes.

Together, the `metablastr` package may enable a new level of data-driven genomics
research by providing the computational tools and data science standards needed
to perform reproducible research at scale.

I happily welcome anyone who wishes to contribute to this project :) Just drop me an email.

### Install `metablastr`

### For Linux Users:

Please install the `libpq-dev` library on you linux machine by typing into the terminal:

```
sudo apt-get install libpq-dev
```

### For all systems install `metablastr` by typing

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("Biostrings", "GenomicFeatures", "GenomicRanges", "Rsamtools", "IRanges", "rtracklayer", "biomaRt"))

# install.packages("devtools")
# install the current version of metablastr on your system
devtools::install_github("HajkD/metablastr", build_vignettes = TRUE, dependencies = TRUE)
```

__Please follow the [Installation Vignette](https://hajkd.github.io/metablastr/articles/installation.html) to install all standalone sequence search tools.__

### Quick start
 
```r
library(metablastr)
# run blastn (nucleotide to nucleotide search) between example query and subject sequences
blast_test <- blast_nucleotide_to_nucleotide(
                 query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
                 subject = system.file('seqs/sbj_nn.fa', package = 'metablastr'),
                 output.path = tempdir(),
                 db.import  = FALSE)
                 
# look at BLAST results
blast_test
```

```
   query_id     subject_id perc_identity num_ident_match… alig_length mismatches
   <chr>        <chr>              <dbl>            <int>       <int>      <int>
 1 333554|PACi… AT1G01010…          84.2              640         760         63
 2 333554|PACi… AT1G01010…          84.0              536         638         90
 3 333554|PACi… AT1G01010…          78.6               44          56         12
 4 470181|PACi… AT1G01020…          94.7              699         738         39
 5 470180|PACi… AT1G01030…          95.2             1025        1077         40
 6 333551|PACi… AT1G01040…          96.0             3627        3779        125
 7 333551|PACi… AT1G01040…          95.5             1860        1948         82
 8 909874|PACi… AT1G01050…          96.6              617         639         22
 9 470177|PACi… AT1G01060…          92.8             1804        1944        110
10 918864|PACi… AT1G01070…          95.3             1046        1098         40
# ... with 13 more rows, and 15 more variables: gap_openings <int>,
#   n_gaps <int>, pos_match <int>, ppos <dbl>, q_start <int>, q_end <int>,
#   q_len <int>, qcov <int>, qcovhsp <int>, s_start <int>, s_end <dbl>,
#   s_len <dbl>, evalue <dbl>, bit_score <chr>, score_raw <int>
```


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/HajkD/metablastr/issues


## Interfaces implemented in `metablastr`:

### Perform BLAST searches 

- `blast_protein_to_protein()`: Perform Protein to Protein BLAST Searches (BLASTP)
- `blast_nucleotide_to_nucleotide()`: Perform Nucleotide to Nucleotide BLAST Searches (BLASTN)
- `blast_nucleotide_to_protein()`: Perform Nucleotide to Protein BLAST Searches (BLASTX)
- `blast_protein_to_nucleotide()`: Perform Protein to Nucleotide BLAST Searches (TBLASTN)
- `blast_best_hit()`: Retrieve only the best BLAST hit for each query
- `blast_best_reciprocal_hit()`: Retrieve only the best reciprocal BLAST hit for each query
- `blast_rpsblast`: Perform Reverse PSI-BLAST searches (rpsblast)
- `read_blast()`: Import BLAST output into R session (in memory) or via `PostgresSQL` database connection.

### BLAST against common NCBI databases 

- `blast_protein_to_nr_database()`: Perform Protein to Protein BLAST Searches against the `NCBI non-redundant database`
- `blast_nt()`: Perform Nucleotide to Nucleotide BLAST Searches against the `NCBI non-redundant database`
- `blast_est()`: Perform Nucleotide to Nucleotide BLAST Searches against the `NCBI expressed sequence tags database`
- `blast_pdb_protein()`:
- `blast_pdb_nucleotide()`:
- `blast_swissprot()`:
- `blast_delta()`:
- `blast_refseq_rna()`:
- `blast_refseq_gene()`:
- `blast_refseq_protein()`:


### BLAST against a set of organisms

- `blast_nucleotide_to_genomes()`: Perfrom BLAST Searches Against a Set of Genomes
- `blast_protein_to_proteomes()`: Perfrom BLAST Searches Against a Set of Proteomes
- `detect_homologs_cds_to_cds()`: Perform CDS to CDS BLAST Searches against a set of CDS files
- `detect_homologs_proteome_to_proteome()`: Perform Proteome to Proteome BLAST Searches against a set of Proteomes
- `extract_hit_seqs_from_genomes()`: Extract sequences of BLAST hits in respective genomes and store it as 'fasta' file(s)
- `extract_random_seqs_from_genome()`: Extract random loci from a genome of interest
- `sample_chromosome_intervals()`: Helper function to sample random intervals of length 'interval_width' from chromosomes

### Analyze BLAST Report

- `filter_blast_`:

### Navigation functions
- `list_outformats()`: List available BLAST output formats

