#FROM dynverse/dynwrappy:v0.1.0
FROM continuumio/miniconda3

ENV SHELL bash

RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get clean

#Install stream package
RUN conda install libgfortran stream -y && conda clean --all -y

COPY definition.yml run.py example.sh /code/

ENTRYPOINT ["/code/run.py"]