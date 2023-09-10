# metablastr

## Seamless Integration of [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) Sequence Searches in R 

### Motivation 

With the rapid expansion of available sequences in biological databases, the landscape of modern life science research is being transformed. Currently, around several hundred thousand genomic sequences from a diverse array of species in the tree of life are freely accessible to the public. The R package [biomartr](https://github.com/ropensci/biomartr) enables automated access and retrieval of this vast data, paving the way to delve into the rich tapestry of sequence diversity, uncovering new insights into evolvability, variation, and the emergence of diseases.

The [biomartr](https://github.com/ropensci/biomartr) package streamlines the retrieval of a massive amount of biological sequence data in a standardized and reproducible manner. Complementing it, the `metablastr` package is tailored to facilitate large-scale sequence searches, also in a standardized and reproducible approach.

In synergy, [biomartr](https://github.com/ropensci/biomartr) and `metablastr` provide researchers with a comprehensive toolset, allowing them to efficiently gather thousands of biological sequences (genomes, proteomes, annotations, etc.) and conduct extensive sequence comparisons using the gold standard sequence search engine [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). This facilitates the extraction of novel patterns highlighting similarities and divergences among vast sets of species.

It's worth noting that the go-to instrument for large-scale sequence searches is BLAST (Basic Local Alignment Search Tool). It is purposefully built to identify regions of sequence similarity between a given query and subject sequences or sequence databases. 

Building on these advancements, we have recently introduced [DIAMOND2](https://www.nature.com/articles/s41592-021-01101-x), a groundbreaking software solution designed to accelerate `BLAST` searches by an factor of up to 10,000x. To offer researchers even more flexibility and integration, we provide [rdiamond](https://github.com/drostlab/rdiamond), a dedicated interface package that allows programmatic handling of [DIAMOND2](https://github.com/bbuchfink/diamond) sequence searches directly through R. This not only streamlines the sequence search process but also ensures that researchers can access and utilize the power of [DIAMOND2](https://github.com/bbuchfink/diamond) within a familiar R environment.

### Short package description  

The `metablastr` package harnesses the power of [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) by providing interface functions between R and the standalone (command line tool) version of this program. 

Together, the `metablastr` package may enable a new level of data-driven genomics
research by providing the computational tools and data science standards needed
to perform reproducible research at scale.

### Citation

__The `metablastr` package is still under development and not formally published yet. However, we did develop parts of `metablastr` for this Methods chapter which you can cite until a `metablastr` specific manuscript is prepared.__

> M Benoit, HG Drost. [__A Predictive Approach to Infer the Activity and Natural Variation of Retrotransposon Families in Plants__](https://link.springer.com/protocol/10.1007%2F978-1-0716-1134-0_1). In: Cho J. (eds) Plant Transposable Elements. Methods in Molecular Biology, vol 2250. Humana, New York, NY (2021). 


### Install `metablastr`

### For Linux Users:

Please install the `libpq-dev` library on you linux machine by typing into the terminal:

```
sudo apt-get install libpq-dev
```

### For all systems install `metablastr` by typing

```r
# install BiocManager if required
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# install package dependencies
BiocManager::install(
  c(
    "Biostrings",
    "GenomicFeatures",
    "GenomicRanges",
    "Rsamtools",
    "IRanges",
    "rtracklayer")
)

# install.packages("devtools")
# install the current version of metablastr on your system
devtools::install_github("drostlab/metablastr", build_vignettes = TRUE, dependencies = TRUE)
```

__Please follow the [Installation Vignette](https://drostlab.github.io/metablastr/articles/installation.html) to install all standalone sequence search tools.__

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

https://github.com/drostlab/metablastr/issues


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

