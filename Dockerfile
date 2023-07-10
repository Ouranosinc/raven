# vim:set ft=dockerfile:
FROM condaforge/mambaforge
ARG DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore
LABEL org.opencontainers.image.authors="smith.trevorj@ouranos.ca"
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.18.2"

WORKDIR /opt/wps
COPY . /opt/wps

# Update Debian system
RUN apt-get update
RUN apt-get install -y build-essential
RUN rm -rf /var/lib/apt/lists/*

# Update conda and mamba
RUN mamba update -n base conda mamba

# Create conda environment
COPY environment.yml .
RUN mamba env create -n raven -f environment.yml
RUN mamba clean --all --yes

# Add the raven conda environment to the path
ENV PATH /opt/conda/envs/raven/bin:$PATH

# Install RavenWPS in editable mode
RUN make install

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099
CMD ["exec raven-wps start -b 0.0.0.0 -c opt/wps/etc/demo.cfg"]

# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
