

:warning: <span style="color:red"> **THIS PROJECT HAS BEEN MOVED TO GITHUB:
[RFLOMICS](https://github.com/RFLOMICS/RFLOMICS)** </span> :warning:

# RFLOMICS: R package and Shiny interface for Integrative analysis of omics data

The acquisition of **multi-omics data** in the context of **complex experimental design** is a widely used practice to identify features (transcripts, proteins, metabolites,...) and decipher the biological processes they are involved in. The investigation of each omics layer is a good first step to explore and extract relevant biological variability. The statistical integration could then be restrained to pertinent omics levels and features. Such analysis of heterogeneous data remains a technical challenge with the needs of expertise methods and parameters to take into account data specificity. Furthermore, applying different statistical methods from several tools is also a technical challenge in term of data management. In this context, we developed RFLOMICS.

**RFLOMICS** is an R package and Shiny interface that offers guided, comprehensive, and reproducible analysis within a framework designed to manage multiple omics datasets. The interface provides remarkable flexibility, enabling users to seamlessly navigate between result exploration and visualization.

<img src="man/figures/workflow.png" align="center" width="600"/>

**Specifications:**

- Performs complete multi-omics project analysis,
- Support multi-factorial experimental design (up to 3 biological factors), 
- Guarantees the relevance of the used methods,
- Reduces the analysis time on the unavoidable steps,
- Ensure the reproducibility of omics analysis,
- Accessible via one and simple user-friendly interface

**Constraints:**

- 3 types of omics data : RNAseq (read count), metabolomics (abundance) and proteomics (abundance)
- less than 10 datasets per omics type.
- 1 to 3 biological factors
- 1 to 2 batch factors and at least 3 replicates per condition
- Complete (mandatory) and balanced (recommended) design 

**Repositories:**

- [forgeMIA](https://forgemia.inra.fr/flomics/rflomics/): official repository
- [github](https://github.com/ijpb-bioinformatics/RFLOMICS/): mirror repository


## Use Rflomics

### SK8 platform
You can access the rflomics app on Sk8 using this url: 

https://rflomics.sk8.inrae.fr/

WARNING: only 2Gb of ram available per app, we  recommand using the SK8 solution 
only if your data have a small number of samples and features. 

### Locally 

#### Install from archive 

Download the archive from 
[forgeMIA](https://forgemia.inra.fr/flomics/rflomics/-/archive/master/rflomics-master.tar.gz)
(tar.gz) or from the 
[github mirror](https://github.com/ijpb-bioinformatics/RFLOMICS/archive/refs/heads/master.zip) (zip)

You can then install the package from the archive using either the install 
utility from Rstudio or running the following command:

``` r
install.packages("rflomics-master.tar.gz", repos = NULL, type = "source")
```

#### Install from forgeMIA or github using remotes 
In a R console, use the following command to install from either the forgeMIA
repository or the github repository:

``` r
remotes::install_gitlab(repo = "flomics/rflomics", host = "https://forgemia.inra.fr/")
```

``` r
remotes::install_github("https://github.com/ijpb-bioinformatics/RFLOMICS")
```


#### Install with git repository 

You can clone the repository from forgeMIA using the following command in a 
shell terminal: 

```  
git clone https://forgemia.inra.fr/flomics/rflomics.git
```

In a R console, set the working directory to the git repository, then install
the package:

``` r
setwd("rflomics/")
remotes::install_local(upgrade="never")
```

### Troubleshooting MOFA2

Rflomics uses the R package [MOFA2](https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html). 
This package depends on a python script and this
can lead to several issues. The first step is to consult the MOFA2 troubleshooting
[FAQ](https://biofam.github.io/MOFA2/troubleshooting.html) on their website.

If none of the proposed solutions are resolving your issues, you can try some
additional steps and verifications. 

1. Install [Python](https://www.python.org/downloads/) or make sure you have it installed.

You can check you have python installed by running `which python` in a shell terminal. 
On windows, replace with `where python`. You should be able to see the path(s) were
Python binaries are installed. You can have multiple python version installed at once.
You have to make sure the one used by R is this correct one. 

2. Install pip command (if not already done). For windows users, you can use [this guide](https://phoenixnap.com/kb/install-pip-windows).
3. Install mofapy 2: in a shell terminal, run `pip install mofapy2`.
4. Set the value of the RETICULATE_PYTHON environment variable to a Python 
binary and restart your session:
* either manually indicate it in your .Rprofile (hidden file) 
* or use the following command in a terminal, in the .Rprofile folder:

```
echo "Sys.setenv(RETICULATE_PYTHON = \"path_to_python_bin\")" >> .Rprofile
```
5. Install MOFA2 using `BiocManager::install("MOFA2")` in R. 

You can check your configuration using `reticulate::py_config()` in R. 
Your Python binary path and mofapy2 path should appear after step 5. 
If it's not the case, something went wrong. 


## Run Rflomics

To run Rflomics, use the runRFLOMICS() function:

``` r
library(RFLOMICS)
runRFLOMICS()
```

## For additional information: [vignettes](https://flomics.pages.mia.inra.fr/rflomics/index.html)

## Licence
Artistic-2.0


## Contact and support
[ijpb-bioinfo-team](mailto:ijpb-bioinfo-team@inrae.fr) (ijpb-bioinfo-team@inrae.fr)

## References
-   [CATI Sysmics](https://sysmics.cati.inrae.fr/),
-   [Ilana L. et al. (2020)](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=32426025&retmode=ref&cmd=prlinks)
-   [Efstathiou G, et al. 2017]()
