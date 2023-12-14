FROM rocker/r-ver:4.3.2

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
    libglpk-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /home

RUN git clone --branch  sk8-test  https://forgemia.inra.fr/flomics/rflomics.git

WORKDIR /home/rflomics

COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

RUN R -e "renv::restore()"

RUN Rscript -e 'remotes::install_local(".",upgrade="never")'

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" >> /usr/local/lib/R/etc/Rprofile.site

RUN addgroup --system rfuser \
    && adduser --system --ingroup rfuser rfuser

#RUN chown -R rfuser:rfuser /home/rfuser
# RUN chown -R rfuser:rfuser .

# USER rfuser

#WORKDIR /home/rfuser

# Install mofapy2
RUN pip install mofapy2==0.7.0


EXPOSE 3838

CMD ["R","-e","reticulate::use_python(\"/usr/bin/python3\", required = NULL);library(RFLOMICS); RFLOMICS::runRFLOMICS()"]



