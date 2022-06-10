# Set the base image to Debian
FROM debian

# File Author / Maintainer
LABEL authors="Romain Daveau <romain.daveau@newbiologix.com>"

# Begin installation
RUN apt-get update && \
    apt-get install -y --no-install-suggests --no-install-recommends perl \
    r-cran-optparse \
    r-cran-rcolorbrewer \
    r-cran-ggplot2 \
    r-cran-reshape2 \
    r-cran-stringr && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# File copy and mode
COPY snptyping-1.pl /usr/local/bin/
COPY snptyping-2.pl /usr/local/bin/
COPY snptyping-2.R /usr/local/bin/
COPY snptyping-3.R /usr/local/bin/
COPY sprt_summary-2.pl /usr/local/bin/
COPY sprt_summary-3.pl /usr/local/bin/

COPY snptyping_rxli-1.pl /usr/local/bin/
COPY snptyping_rxli-2.pl /usr/local/bin/
COPY snptyping_rxli-2.R /usr/local/bin/
COPY snptyping_rxli-3.R /usr/local/bin/
COPY sprt_summary_rxli-2.pl /usr/local/bin/
COPY sprt_summary_rxli-3.pl /usr/local/bin/

RUN chmod 755 /usr/local/bin/snptyping*
RUN chmod 755 /usr/local/bin/sprt*

# Change default user
RUN useradd -ms /bin/bash docker
USER docker
