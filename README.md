# metablastr
### Perform Local BLAST Searches with R
The Basic Local Alignment Search Tool (BLAST) finds regions of sequence similarity between a query and a subject sequence or sequence database.

The `metablastr` package aims to provide R interface functions for the standalone (command line tool) version
of BLAST. This way users can perform local BLAST searches with R and import the results back
to the R session.


```r
# install.packages("devtools")

# install the current version of blastr on your system
library(devtools)
install_github("HajkD/metablastr", build_vignettes = TRUE, dependencies = TRUE)
```


## Interfaces implemented in `metablastr`:

### Perform BLAST searches 
- `blast_protein()`: Perform Protein to Protein BLAST Searches


### Navigation functions
- `list_outformats()`: List available BLAST output formats

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/HajkD/metablastr/issues
