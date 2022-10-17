FROM rocker/tidyverse:4.1.2

USER root

WORKDIR /home/rstudio

# Install libs
RUN apt-get update && \
	apt-get -y --no-install-recommends install --fix-missing \
        libglpk-dev \
        libxt-dev \
        openjdk-11-jdk \
        ca-certificates-java \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

RUN wget https://hiv.bio.ed.ac.uk/clusterPicker/ClusterPicker_1.2.3.jar

# Install R packages
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('tidytree','ape','parallel'), repos='http://cran.rstudio.com')"

RUN R -e "BiocManager::install(c('treeio', 'ggtree', 'ggtreeExtra'))"
