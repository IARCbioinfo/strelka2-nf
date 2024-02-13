################## BASE IMAGE #####################
FROM continuumio/miniconda3:23.10.0-1

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="23.10.0-1"
LABEL software="strelka2-nf"
LABEL software.version="1.3"
LABEL about.summary="Container image containing all requirements for strelka2-nf"
LABEL about.home="http://github.com/IARCbioinfo/strelka2-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/strelka2-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/strelka2-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@iarc.fr**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n strelka2-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/strelka2-nf/bin:/opt/conda/envs/strelka2-nf/share/strelka-2.9.10-1/bin:$PATH
