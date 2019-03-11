# MoDentify
Phenotype-driven module identification based on multilevel metabolomics networks

R package

**!!! Broken Cytoscape connection: !!!**

We are aware that the connection to Cytoscape does not work anymore in the current version of MoDentify. We are working on this issue and hope to have a fix soon.

In the meantime, the connection should work with Bioconductor Version 3.6 or older.

- Jan Krumsiek, February 13, 2019.

**Citation:**

Do KD, Rasp D, Kastenmueller G, Suhre K, Krumsiek J.  
MoDentify: phenotype-driven module identification in metabolomics networks at different resolutions.  
Bioinformatics, 2018  
https://doi.org/10.1093/bioinformatics/bty650

**Installation:**

```r
# Installation of the development version from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("David-J-R/MoDentify", build_opts = c())
```

(takes about a minute)

**Show vignette with example code:**

`browseVignettes("MoDentify")`
