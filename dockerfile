# Dockerfile for RUMP

FROM nfcore/base

LABEL maintainer="xinsongdu@ufl.edu"

USER root

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
    vim \
    procps \
    libfreetype6 \
    libcairo2-dev \
    libexpat1-dev \
    libgmp3-dev \
    liblapack-dev \
    libnetcdf-dev \
    libopenbabel-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libgsl0-dev \
    libmpfr-dev \
    pkg-config \
    fftw3-dev \
    libgtk2.0-dev \
    libtiff5-dev \
    libnetcdf-dev \
    libmpfr-dev \
    libnetcdf-dev \
    liblapack-dev \
    cmake \
    python\
    python-dev\
    software-properties-common\
    python-pip\
    python3-pip\
    python-tk\
    python3-tk\
    libnetcdf-dev libpng-dev libbz2-dev liblzma-dev libpcre3-dev libicu-dev

# Install python3-based necessary dependencies for UMPIRE
# Install the conda environment
COPY environment.yml /
RUN conda env update -n rump -f/environment.yml && conda clean -a

RUN pip install mummichog

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/rump/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name rump > /rump.yml

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

# build dirs for UFRC
RUN mkdir /ufrc /orange /bio /rlts 
RUN mkdir -p /scratch/local
RUN mkdir app

# define work dir
WORKDIR /app
COPY accessibility.properties /app

# Fix a bug for java
# RUN mv accessibility.properties /etc/java-8-openjdk/

# install R packages
COPY r_package_install.R /app
RUN Rscript /app/r_package_install.R

# Adds scripts inside the container
RUN mkdir /opt/rump
COPY ./rump/ /opt/rump

# Make files within rump executable
RUN find /opt/rump/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/rump/ -type f -iname "*.R"   -exec chmod +x {} \;

# # Add MZmine in container
# RUN wget https://github.com/mzmine/mzmine2/releases/download/v2.53/MZmine-2.53-Linux.zip && \
#     unzip MZmine-2.53-Linux.zip -d / && \
#     rm MZmine-2.53-Linux.zip && \
#     mv /MZmine-2.53-Linux /MZmine

# Add rump folder with scripts to PATH
ENV PATH /opt/rump:$PATH
