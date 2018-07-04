# vim:set ft=dockerfile:
FROM continuumio/miniconda3
MAINTAINER https://github.com/huard/raven
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.1.0"

# Update Debian system
RUN apt-get update && apt-get install -y \
 build-essential \
&& rm -rf /var/lib/apt/lists/*

# Update conda
RUN conda update -n base conda

# Copy WPS project
COPY . /opt/wps

WORKDIR /opt/wps

# Create conda environment
RUN conda env create -n wps -f environment.yml

# Install WPS
RUN ["/bin/bash", "-c", "source activate wps && python setup.py develop"]

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["source activate wps && exec emu start -b 0.0.0.0 -config /opt/wps/etc/demo.cfg"]

# docker build -t huard/raven .
# docker run -p 9099:9099 huard/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
