FROM rocker/shiny:4.5.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    zip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/shiny-server
COPY install_dependencies.R /tmp/install_dependencies.R
RUN Rscript /tmp/install_dependencies.R

COPY . /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server

USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
