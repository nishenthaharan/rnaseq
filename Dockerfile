FROM rocker/shiny:4.5.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
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

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
  CMD curl --fail --silent http://localhost:3838/ > /dev/null || exit 1

CMD ["/usr/bin/shiny-server"]
