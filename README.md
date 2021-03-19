# MoDentify
Phenotype-driven module identification based on multilevel metabolomics networks

R package

**Version history**
- 3/11/19 - Updated package to work with RCy3 package for Cytoscape connectivity
- 1/18/18 - Initial commit

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
devtools::install_github("krumsieklab/MoDentify", build_vignettes = T)
```

(takes about a minute, since it compiles the Vignette)

**Show vignette with example code:**

`browseVignettes("MoDentify")`
