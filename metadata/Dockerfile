FROM rocker/tidyverse:4.1.2

USER root

WORKDIR /home/rstudio

# Install libs
RUN apt-get update && \
	apt-get -y --no-install-recommends install --fix-missing \
        libglpk-dev \
        libxt-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidytree', repos='http://cran.rstudio.com')"

RUN R -e "BiocManager::install('treeio, ggtree, ggtreeExtra')"
