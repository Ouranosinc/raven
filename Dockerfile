# vim:set ft=dockerfile:
FROM condaforge/miniforge3
ARG DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore
LABEL org.opencontainers.image.authors="https://github.com/Ouranosinc/raven"
LABEL org.opencontainers.image.title="Raven WPS"
LABEL org.opencontainers.image.vendor="Birdhouse"
LABEL org.opencontainers.image.version="0.19.0"

# Set the working directory to /code
WORKDIR /code

# Create conda environment
COPY environment.yml .
RUN mamba env create -n raven -f environment.yml && \
    mamba install -n raven --yes gunicorn && \
    mamba clean --all --yes

# Add the raven conda environment to the path
ENV PATH="/opt/conda/envs/raven/bin:$PATH"

# For pyproj, to avoid error "PROJ: proj_create_from_database: Open of /opt/conda/envs/raven/share/proj failed"
ENV PROJ_DATA="/opt/conda/envs/raven/share/proj"

# Copy WPS project
COPY . /code

# Install WPS project
RUN conda run -n raven pip install --no-cache-dir . --no-deps

# Start WPS service on port 9099
EXPOSE 9099

# Specify a non-root user to run the application
RUN useradd --create-home --shell /bin/bash --uid 1001 nonroot && chown -R nonroot:nonroot /code /home/nonroot /opt/conda/envs/raven
USER nonroot

CMD ["gunicorn", "--bind=0.0.0.0:9099", "-t 60", "raven.wsgi:application"]
# docker build -t pavics/raven .
# docker run -p 9099:9099 pavics/raven
# http://localhost:9099/wps?request=GetCapabilities&service=WPS
# http://localhost:9099/wps?request=DescribeProcess&service=WPS&identifier=all&version=1.0.0
