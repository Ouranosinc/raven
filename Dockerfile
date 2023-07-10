# vim:set ft=dockerfile:
FROM condaforge/mambaforge
ARG DEBIAN_FRONTEND=noninteractive
MAINTAINER Trevor James Smith <smith.trevorj@ouranos.ca>
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.18.2"

# Update Debian system
RUN apt-get update && apt-get install -y build-essential unzip libnetcdf-dev curl && rm -rf /var/lib/apt/lists/*

# Update conda and mamba
RUN mamba update -n base conda mamba

# Create conda environment
COPY environment.yml /opt/wps/
RUN mamba create --yes -n wps python=3.9 && mamba env update -n wps -f /opt/wps/environment.yml

# Copy WPS project
COPY . /opt/wps

# Switch to WPS project directory
WORKDIR /opt/wps

# Install RavenWPS in editable mode
RUN ["/bin/bash", "-c", "source activate wps && make install"]

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099
CMD ["/bin/bash", "-c", "source activate wps && exec raven-wps start -b 0.0.0.0 -c /opt/wps/etc/demo.cfg"]

# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
