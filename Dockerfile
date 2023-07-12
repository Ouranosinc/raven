# vim:set ft=dockerfile:
FROM condaforge/mambaforge
ARG DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore
LABEL org.opencontainers.image.authors="smith.trevorj@ouranos.ca"
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.18.2"

WORKDIR /code
COPY . /code

# Create conda environment
RUN mamba env create -n raven -f environment.yml
# Remove conda cache
RUN mamba clean --all --yes

# Add the raven conda environment to the path
ENV PATH /opt/conda/envs/raven/bin:$PATH

# Install RavenWPS
RUN pip install . --no-deps

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099

# CMD["gunicorn", "--bind=0.0.0.0:5000", "-t 60", "finch.wsgi:application"]
CMD ["exec raven-wps start -b 0.0.0.0 -c /code/etc/demo.cfg"]

# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
