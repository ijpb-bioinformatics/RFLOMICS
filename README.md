# RFLOMICS : R package and Shiny interface for Integrative analysis of omics data

## Presentation

The acquisition of **multi-omics data** in the context of **complex
experimental design** is a widely used practice to identify genes,
proteins, metabolites and decipher the biological processes they are
involved in. The investigation of each omics layer is a good first step
to explore and extract relevant biological variability. The statistical
integration could then be restrained to pertinent omics levels and
features. Such analysis of heterogeneous data remains a technical
challenge with the needs of expert methods and parameters to take into
account data specificity. Furthermore, applying different statistical
methods from several tools is also a technical challenge in term of data
management. In this context, we developed RFLOMICS: **R package coupled
with a shiny application** dedicated to the management and analysis of
multiple omics-datasets in the statistical framework of **vertical
integration** of observations (i.e. analysis of omics data across
experiments on the same individuals) see the figure below.


RFLOMICS currently supports up to three types of omics: RNAseq,
proteomics, and metabolomics.

<img src="www/workflow.png" align="center" width="600"/>

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
- 1 - 3 biological factors
- 1 - 2 batch factors and at least 3 replicates
- Complete (mandatory) and balanced (recommended) design 


## Use RFLOMICS locally

### Install

Download from <https://forgemia.inra.fr/flomics/rflomics/-/tree/v0.1>

``` r
install.packages("rflomics.tar.gz", repos = NULL, type = "source")
```

Or

Clone from forgemia repository

```         
git clone -branch  v0.2  https://forgemia.inra.fr/flomics/rflomics.git
```

``` r
library(remotes)

setwd("rflomics/")

remotes::install_local(upgrade="never")
```

### Run rflomics

``` r
library(RFLOMICS)

RFLOMICS::runRFLOMICS()
```

### Install via Docker

````{=html}
<!--## Install RFLOMICS via Docker

* install Docker

https://docs.docker.com/engine/install/

* Get the Dockerfile

```
https://forgemia.inra.fr/flomics/rflomics/-/blob/develop.0.1/Dockerfile
```

* Build the Docker Image

```
docker build --file=Dockerfile --tag=rflomics .
```

* Run Docker

```
docker run -it -p 3838:3838 -v ${HOME}:/root --name='rflomics' --cpus 4 rflomics
```

* Open a web navigator and paste this url:

```
http://0.0.0.0:3838
```
-->
````

### [Vignettes](https://flomics.pages.mia.inra.fr/rflomics/index.html)

## Licence
blabla

## Contact and support

[ijpb-bioinfo-team](mailto:ijpb-bioinfo-team@inrae.fr) (ijpb-bioinfo-team@inrae.fr)

## References

-   [CATI Sysmics](https://sysmics.cati.inrae.fr/),
-   [Ilana L. et al. (2020)](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=32426025&retmode=ref&cmd=prlinks)
-   [Efstathiou G, et al. 2017]()
