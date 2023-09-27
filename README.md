# RFLOMICS : R package and Shiny interface for Integrative analysis of omics data

## Presentation

<p style="text-align: justify;">
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
</p>

RFLOMICS currently supports up to three types of omics: RNAseq,
proteomics, and metabolomics.

<img src="man/figures/Rflomics_features.png" align="center" width="600"/>

**Aims** - Guarantee the relevance of the used methods and parameters
(RNAseq workflow:
[DicoExpress](https://plantmethods.biomedcentral.com/articles/10.1186/s13007-020-00611-7),
[CATI Sysmics](https://sysmics.cati.inrae.fr/)) - Ensure the
reproducibility of analysis

**Features** - It can deal with **multi-factorial experiments** (up to 3
biological factors). - It can allows the remote computing for time/cpu
consuming tasks **(clustermq package)**

## Locally installation (single omics analysis)

Download from <https://forgemia.inra.fr/flomics/rflomics/-/tree/v0.1>

``` r
install.packages("rflomics.tar.gz", repos = NULL, type = "source")
```

Or

Clone from forgemia repository

```         
git clone -branch  v0.1  https://forgemia.inra.fr/flomics/rflomics.git
```

``` r
library(remotes)

setwd("rflomics/")

remotes::install_local(upgrade="never")
```

## Run rflomics

``` r
library(RFLOMICS)

RFLOMICS::runRFLOMICS()
```

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

## Contact and support

[ijpb-bioinfo-team](mailto:ijpb-bioinfo-team@inrae.fr)

## References

-   [CATI Sysmics](https://sysmics.cati.inrae.fr/),
-   [Ilana L. et al. (2020),
    DiCoExpress](http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=32426025&retmode=ref&cmd=prlinks)
-   [MultiAssayExperiment
    package](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)
-   [edgeR
    package](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
-   [limma
    package](https://bioconductor.org/packages/release/bioc/html/limma.html)
-   [coseq
    package](https://bioconductor.org/packages/release/bioc/html/coseq.html)
-   [clustermq
    package](https://cran.r-project.org/web/packages/clustermq/index.html)
