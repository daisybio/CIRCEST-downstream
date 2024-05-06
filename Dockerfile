FROM rocker/r2u:jammy

WORKDIR /app
COPY setup.R .
RUN Rscript setup.R

COPY run.R src data .

CMD ["Rscript", "run.R"]