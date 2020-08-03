# Dockerfile for RUMP

FROM rocker/rstudio:3.6.3

LABEL maintainer="xinsongdu@ufl.edu"

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
    vim \
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
    default-jdk \
    python\
    python-dev\
    software-properties-common\
    python-pip\
    python3-pip\
    python-tk\
    python3-tk\
    libnetcdf-dev libpng-dev libbz2-dev liblzma-dev libpcre3-dev libicu-dev

# Install python3-based necessary dependencies for UMPIRE
RUN pip3 install --upgrade 'setuptools==45.2.0'
RUN pip3 install 'numpy==1.18.1' 'scipy==1.4.1' 'pandas==0.25.3' 'matplotlib<3.0.0,>=2.1.1' 'plotly==4.5.0' 'seaborn==0.9.1' 'scikit-learn==0.22.1' matplotlib_venn 'multiqc==1.8' 'statsmodels==0.11.0' 'fastcluster==1.1.26' 'pylint==2.4.4'
RUN echo "alias python=python3" >> ~/.bash_profile

ENV NETCDF_INCLUDE=/usr/include

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

# Install mummichog
RUN pip install --upgrade 'setuptools==44.0.0'
RUN pip install 'mummichog==2.2.0'

# Adds scripts inside the container
RUN mkdir /opt/rump
COPY ./rump/ /opt/rump

# Make files within rump executable
RUN find /opt/rump/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/rump/ -type f -iname "*.R"   -exec chmod +x {} \;

# Add MZmine in container
RUN wget https://github.com/mzmine/mzmine2/releases/download/v2.53/MZmine-2.53-Linux.zip && \
    unzip MZmine-2.53-Linux.zip -d / && \
    rm MZmine-2.53-Linux.zip && \
    mv /MZmine-2.53-Linux /MZmine

# Add rump folder with scripts to PATH
ENV PATH /opt/rump:$PATH
