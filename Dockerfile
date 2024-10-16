# vim:set ft=dockerfile:
FROM condaforge/mambaforge
ARG DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore
LABEL org.opencontainers.image.authors="https://github.com/Ouranosinc/raven"
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.18.2"

# Set the working directory to /code
WORKDIR /code

# Create conda environment
COPY environment.yml .
RUN mamba env create -n raven -f environment.yml && mamba install -n raven gunicorn  && mamba clean --all --yes

# Add the raven conda environment to the path
ENV PATH /opt/conda/envs/raven/bin:$PATH

# Copy raven source code
COPY . /code

# Install raven
RUN pip install . --no-deps

# Start WPS service on port 9099 on 0.0.0.0
EXPOSE 9099

CMD ["gunicorn", "--bind=0.0.0.0:9099", "raven.wsgi:application"]
#CMD ["exec raven-wps start -b '0.0.0.0' -c etc/demo.cfg"]

# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
