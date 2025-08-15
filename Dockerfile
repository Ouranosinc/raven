# vim:set ft=dockerfile:
FROM condaforge/mambaforge
ARG DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore
LABEL org.opencontainers.image.authors="https://github.com/Ouranosinc/raven"
LABEL Description="Raven WPS" Vendor="Birdhouse" Version="0.19.0"

# Set the working directory to /code
WORKDIR /code

# Create conda environment
COPY environment.yml .
RUN mamba env create -n raven -f environment.yml && mamba install -n raven gunicorn && mamba clean --all --yes

# Add the raven conda environment to the path
ENV PATH=/opt/conda/envs/raven/bin:$PATH

# Copy WPS project
COPY . /code

# Install WPS project
RUN pip install . --no-deps

# Set environment variables for Raven WPS
ENV RAVEN_BIND_HOST=0.0.0.0
ENV RAVEN_BIND_PORT=9099

# Start WPS service on port 9099
EXPOSE 9099

CMD ["sh", "-c", "gunicorn --bind=${RAVEN_BIND_HOST}:${RAVEN_BIND_PORT} -t 60 raven.wsgi:application"]
# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
