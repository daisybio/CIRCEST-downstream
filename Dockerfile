FROM rocker/r2u:jammy

WORKDIR /app
COPY setup.R .
RUN Rscript setup.R

COPY run.R .
COPY data data
COPY src src

CMD ["Rscript", "run.R"]