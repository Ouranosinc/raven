# vim:set ft=dockerfile:
FROM continuumio/miniconda3
MAINTAINER https://github.com/huard/raven
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.15.1"

# Update Debian system
RUN apt-get update && apt-get install -y \
 build-essential unzip libnetcdf-dev curl \
&& rm -rf /var/lib/apt/lists/*

# Update conda
RUN conda update -n base conda

# Create conda environment
COPY environment.yml /opt/wps/
RUN conda create --yes -n wps python=3.7 && conda env update -n wps -f /opt/wps/environment.yml

# Copy WPS project
COPY . /opt/wps

WORKDIR /opt/wps

# (1) Install RavenPy: note that this also installs the Raven and Ostrich binaries in the wps conda env's bin
# (2) Install RavenWPS in editable mode
# Have to uninstall the ravenpy installed by conda so the re-install with
# binaries work.  Same problem as in the Makefile.
RUN ["/bin/bash", "-c", "source activate wps && pip install -e ."]

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099
CMD ["/bin/bash", "-c", "source activate wps && exec raven-wps start -b 0.0.0.0 -c /opt/wps/etc/demo.cfg"]

# docker build -t huard/raven .
# docker run -p 9099:9099 huard/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
