FROM rocker/shiny-verse:4.1.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libgmp-dev \
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

RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'

RUN Rscript -e 'BiocManager::install("S4Vectors",version="3.14")'
RUN Rscript -e 'BiocManager::install("edgeR",version="3.14")'
RUN Rscript -e 'BiocManager::install("limma",version="3.14")'
RUN Rscript -e 'BiocManager::install("MultiAssayExperiment",version="3.14")'
RUN Rscript -e 'BiocManager::install("SummarizedExperiment",version="3.14")'
RUN Rscript -e 'BiocManager::install("pheatmap",version="3.14")'
RUN Rscript -e 'BiocManager::install("coseq",version="3.14")'
RUN Rscript -e 'BiocManager::install("BiocParallel",version="3.14")'
RUN Rscript -e 'BiocManager::install("EnhancedVolcano",version = "3.14")'
RUN Rscript -e 'BiocManager::install("ComplexHeatmap",upgrade="never", version = "3.14")'
RUN Rscript -e 'BiocManager::install("MOFA2",upgrade="never", version = "3.14")'


RUN Rscript -e 'remotes::install_version("shinyBS",upgrade="never", version = "0.61")'
RUN Rscript -e 'remotes::install_version("reactable",upgrade="never", version = "0.2.3")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.6.2")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("shinyFiles",upgrade="never", version = "0.9.1")'

RUN Rscript -e 'remotes::install_version("venn",upgrade="never", version = "1.10")'
RUN Rscript -e 'remotes::install_version("gridExtra",upgrade="never", version = "2.3")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.0")'
RUN Rscript -e 'remotes::install_version("UpSetR",upgrade="never", version = "1.4.0")'
RUN Rscript -e 'remotes::install_version("ggpubr",upgrade="never", version = "0.4.0")'


RUN Rscript -e 'remotes::install_version("reshape2",upgrade="never", version = "1.4.4")'
RUN Rscript -e 'remotes::install_version("kableExtra",upgrade="never", version = "1.3.4")'

RUN Rscript -e 'remotes::install_version("clustermq",upgrade="never", version = "0.8.95.2")'
RUN Rscript -e 'remotes::install_version("devtools",upgrade="never", version = "2.4.3")'

RUN Rscript -e 'remotes::install_version("FactoMineR",upgrade="never", version = "2.4")'
RUN Rscript -e 'remotes::install_version("statmod",version="1.4.36")'

WORKDIR /home

RUN git clone --branch  develop  https://forgemia.inra.fr/flomics/rflomics.git

WORKDIR /home/rflomics

RUN Rscript -e 'remotes::install_local(upgrade="never")'

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" >> /usr/local/lib/R/etc/Rprofile.site


RUN addgroup --system rfuser \
    && adduser --system --ingroup rfuser rfuser

RUN chown -R rfuser:rfuser /home/rfuser

USER rfuser

WORKDIR /home/rfuser


EXPOSE 3838

CMD ["R","-e","library(RFLOMICS); RFLOMICS::runExample()"]

