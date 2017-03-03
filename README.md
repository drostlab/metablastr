# blastr
### Performing local BLAST searches with R
The Basic Local Alignment Search Tool (BLAST) finds regions of sequence similarity between a query and a subject sequence or sequence database.

The `blastr` package aims to provide R interface functions for the standalone (command line tool) version
of BLAST. This way users can perform local BLAST searches with R and import the results back
to the R session.


```r
# install.packages("devtools")

# install the current version of blastr on your system
library(devtools)
install_github("HajkD/blastr", build_vignettes = TRUE, dependencies = TRUE)
```


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/HajkD/blastr/issues
