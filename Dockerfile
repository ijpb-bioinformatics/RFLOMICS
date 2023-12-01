FROM rocker/r-ver:4.3.1

LABEL maintainer="ijpb-bioinfo-team@inrae.fr"
LABEL version="0.1"
LABEL description="Dockerfile to construct an iamge of the rflomics application"


RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libgmp-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \ 
    libtiff5-dev \
    libudunits2-dev \
    libxml2-dev \
    libzmq3-dev \
    libnlopt-dev \
    libcairo2-dev \
    libxt-dev \
    libssh2-1-dev \
    openssh-server\
    git-all \
    cmake \
    python3 \
    python3-setuptools \
    python3-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists/*

RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/', Bioc = 'http://www.bioconductor.org/packages/release/bioc'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site


RUN R -e 'install.packages(c("remotes", "BiocManager"))'

RUN R -e 'remotes::install_github("r-lib/remotes", ref = "97bbf81")'

RUN R -e 'BiocManager::install("AnnotationDbi",upgrade="never",version="1.62.2")'
RUN R -e 'BiocManager::install("BiocParallel",upgrade="never",version="1.34.2")'
RUN R -e 'BiocManager::install("clusterProfiler",upgrade="never",version="4.8.3")'
RUN R -e 'BiocManager::install("ComplexHeatmap",upgrade="never",version="2.16.0")'
RUN R -e 'BiocManager::install("coseq",upgrade="never",version="1.24.0")'
RUN R -e 'BiocManager::install("edgeR",upgrade="never",version="edgeR_3.42.4")'
RUN R -e 'BiocManager::install("EnhancedVolcano",upgrade="never",version="1.18.0")'
RUN R -e 'BiocManager::install("enrichplot",upgrade="never",version="1.20.0")'
RUN R -e 'BiocManager::install("limma",upgrade="never",version="3.56.2")'
RUN R -e 'BiocManager::install("mixOmics",upgrade="never",version="6.24.0")'
RUN R -e 'BiocManager::install("MOFA2",upgrade="never",version="1.10.0")'
RUN R -e 'BiocManager::install("MultiAssayExperiment",upgrade="never",version="1.26.0")'
RUN R -e 'BiocManager::install("pathview",upgrade="never",version="1.40.0")'
RUN R -e 'BiocManager::install("pheatmap",upgrade="never",version="1.0.12")'
RUN R -e 'BiocManager::install("S4Arrays",upgrade="never",version="1.0.6")'
RUN R -e 'BiocManager::install("S4Vectors",upgrade="never",version="0.38.2")'
RUN R -e 'BiocManager::install("SummarizedExperiment",upgrade="never",version="1.30.2")'

RUN Rscript -e 'remotes::install_version("clustermq", upgrade="never",version="0.9.0")'
RUN Rscript -e 'remotes::install_version("colourpicker",upgrade="never",version="1.3.0")'
RUN Rscript -e 'remotes::install_version("data.table",upgrade="never",version="1.14.8")'
RUN Rscript -e 'remotes::install_version("DT", upgrade="never",version="0.30")'
RUN Rscript -e 'remotes::install_version("FactoMineR", upgrade="never",version="2.9")'
RUN Rscript -e 'remotes::install_version("gridExtra",upgrade="never",version="2.3")'
RUN Rscript -e 'remotes::install_version("ggpubr", upgrade="never",version="0.6.0")'
RUN Rscript -e 'remotes::install_version("htmltools", upgrade="never",version="0.5.7")'
RUN Rscript -e 'remotes::install_version("kableExtra", upgrade="never",version="1.3.4")'
RUN Rscript -e 'remotes::install_version("qgraph",upgrade="never",version="1.9.8")'
RUN Rscript -e 'remotes::install_version("RColorBrewer",upgrade="never",version="1.1-3")'
RUN Rscript -e 'remotes::install_version("reactable",upgrade="never",version="0.4.4")'
RUN Rscript -e 'remotes::install_version("reshape2", upgrade="never",version="1.4.4")'
RUN Rscript -e 'remotes::install_version("reticulate", upgrade="never",version="1.34.0")'
RUN Rscript -e 'remotes::install_version("rmarkdown", upgrade="never",version="2.25")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version="1.7.5.1")'
RUN Rscript -e 'remotes::install_version("shinyBS",upgrade="never",version="0.61.1")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version="0.7.2")'
RUN Rscript -e 'remotes::install_version("shinyFiles",upgrade="never",version="0.9.3")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never",version="0.8.0")'
RUN Rscript -e 'remotes::install_version("statmod", upgrade="never",version="1.5.0")'
RUN Rscript -e 'remotes::install_version("tidyverse",upgrade="never",version="0.4.5")'
RUN Rscript -e 'remotes::install_version("UpSetR", upgrade="never",version="1.4.0")'
RUN Rscript -e 'remotes::install_version("vroom", upgrade="never",version="1.6.4")'


WORKDIR /home

RUN git clone --branch  develop_0.2  https://forgemia.inra.fr/flomics/rflomics.git

WORKDIR /home/rflomics

RUN Rscript -e 'remotes::install_local(".",upgrade="never")'

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" >> /usr/local/lib/R/etc/Rprofile.site

RUN addgroup --system rfuser \
    && adduser --system --ingroup rfuser rfuser

RUN chown -R rfuser:rfuser /home/rfuser

USER rfuser

WORKDIR /home/rfuser


EXPOSE 3838

CMD ["R","-e","library(RFLOMICS); RFLOMICS::runRFLOMICS()"]



