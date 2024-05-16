################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.12.0

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.8.2"
LABEL software="MSA"
LABEL software.version="2.1"
LABEL about.summary="Container image containing all requirements for MSA nextflow pipeline"
LABEL about.home="https://gitlab.com/s.senkin/MSA/"
LABEL about.documentation="https://gitlab.com/s.senkin/MSA/README.md"
LABEL about.license_file="https://gitlab.com/s.senkin/MSA/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **s.senkin** <**s.senkin@gmail.com**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n MSA -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/MSA/bin:$PATH
