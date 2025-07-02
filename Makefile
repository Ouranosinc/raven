# Configuration
# Determine this makefile's path.
# Be sure to place this BEFORE `include` directives, if any.
THIS_FILE := $(lastword $(MAKEFILE_LIST))
APP_ROOT := $(abspath $(lastword $(MAKEFILE_LIST))/..)
APP_NAME := raven-wps

OS := $(shell uname)

WPS_PORT := 9099
WPS_URL := http://0.0.0.0:$(WPS_PORT)

# If WPS_URL is overridden, this should also be overridden to match.
WPS_OUTPUT_URL := http://localhost:$(WPS_PORT)/outputs

# This will only work on Linux (not macOS/homebrew GDAL)
GDAL_VERSION := $(shell gdal-config --version)

# Used in target refresh-notebooks to make it looks like the notebooks have
# been refreshed from the production server below instead of from the local dev
# instance so the notebooks can also be used as tutorial notebooks.
OUTPUT_URL = https://pavics.ouranos.ca/wpsoutputs/raven

SANITIZE_FILE := https://github.com/Ouranosinc/PAVICS-e2e-workflow-tests/raw/master/notebooks/output-sanitize.cfg

ANACONDA_HOME := $(shell conda info --base 2> /dev/null)

ifeq "$(ANACONDA_HOME)" ""
ANACONDA_HOME := $(HOME)/miniconda3
endif

CONDA := $(shell command -v conda 2> /dev/null)
CONDA_ENV ?= $(APP_NAME)
PYTHON_VERSION = 3.9

# Choose Anaconda installer depending on your OS
ANACONDA_URL = https://repo.anaconda.com/miniconda
UNAME_S := $(shell uname -s)
DOWNLOAD_CACHE = /tmp/

# Additional servers used by notebooks
FINCH_WPS_URL ?= https://pavics.ouranos.ca/twitcher/ows/proxy/finch/wps

# To run tests on local servers, use
# make FINCH_WPS_URL=http://localhost:5000 test-notebooks

FN ?= "Miniconda3-latest-Linux-x86_64.sh"
GDAL_VERSION = "$(shell gdal-config --version)"
ifeq ($(OS),"Darwin")
	FN = "Miniconda3-latest-MacOSX-x86_64.sh"
	GDAL_VERSION = "$(shell gdalinfo --version | awk '{print $2}' | sed 's/.$//')"
endif
ifeq ($(OS),"Windows_NT")
# UNTESTED
	FN = Miniconda3-latest-Windows-x86_64.sh
	GDAL_VERSION = "$(shell gdalinfo --version)"
endif

# end of configuration

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
	else:
		match = re.match(r'^## (.*)$$', line)
		if match:
			help = match.groups()[0]
			print("\n%s" % (help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

.DEFAULT_GOAL := help

help: ## print this help message. (Default)
	@echo "Please use 'make <target>' where <target> is one of:"
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

## Build targets:

install: ## install raven application
ifdef GDAL_VERSION
	@echo "Installing application with GIS..."
	@-bash -c "pip install --no-cache-dir gdal==$(GDAL_VERSION)"
endif
	@-bash -c "pip install -e ."
	@echo "\nStart service with \`make start\` and stop with \`make stop\`."

develop: ## install raven application with development libraries
	@echo "Installing development requirements for tests and docs ..."
	@-bash -c 'pip install -e ".[dev]"'

start: ## start raven service as daemon (background process)
	@echo "Starting application ..."
	@-bash -c "$(APP_NAME) start -d"

stop: ## stop raven service
	@echo "Stopping application ..."
	@-bash -c "$(APP_NAME) stop"

restart: stop start  ## restart raven service
	@echo "Restarting application ..."

status: ## show status of raven service
	@echo "Showing status ..."
	@-bash -c "$(APP_NAME) status"

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	@echo "Removing build artifacts ..."
	@-rm -fr build/
	@-rm -fr dist/
	@-rm -fr .eggs/
	@-find . -name '*.egg-info' -exec rm -fr {} +
	@-find . -name '*.egg' -exec rm -f {} +
	@-find . -name '*.log' -exec rm -fr {} +
	@-find . -name '*.sqlite' -exec rm -fr {} +

clean-pyc: ## remove Python file artifacts
	@echo "Removing Python file artifacts ..."
	@-find . -name '*.pyc' -exec rm -f {} +
	@-find . -name '*.pyo' -exec rm -f {} +
	@-find . -name '*~' -exec rm -f {} +
	@-find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	@echo "Removing test artifacts ..."
	@-rm -fr .pytest_cache

clean-dist: clean  ## remove git ignored files and directories
	@echo "Running 'git clean' ..."
	@git diff --quiet HEAD || echo "There are uncommitted changes! Aborting 'git clean' ..."
	## do not use git clean -e/--exclude here, add them to .gitignore instead
	@-git clean -dfx

clean-docs: ## remove documentation artifacts
	@echo "Removing documentation artifacts ..."
	@-rm -f docs/raven.rst
	@-rm -f docs/modules.rst
	$(MAKE) -C docs clean

lint/flake8: ## check style with flake8
	python -m ruff check src/raven tests
	python -m flake8 --config=.flake8 src/raven tests
	python -m numpydoc lint src/raven/*/*.py
	yamllint --config-file=.yamllint.yaml src/raven

lint/black: ## check style with black
	python -m black --check src/raven tests
	python -m blackdoc --check src/raven docs

lint: lint/flake8 lint/black ## check style

## Testing targets:

test: ## run tests quickly with the default Python
	@echo "Running tests (skip slow and online tests) ..."
	@bash -c "pytest -v -m 'not slow and not online and not very_slow' tests/"

test-advanced: ## run tests with slow and online tests
	@echo "Running advanced tests (skip very_slow tests) ..."
	@bash -c "pytest -v -m 'not very_slow' tests/"

test-all: ## run all tests
	@echo "Running all tests (including slow and online tests) ..."
	@bash -c "pytest -v tests/"

test_pdb: ## run tests quickly with the default Python and drop into pdb
	@echo "Running tests (skip slow and online tests) with --pdb ..."
	@bash -c "pytest -v -m 'not slow and not online' --pdb"

test-tox: ## run tests on every available Python version with tox
	@bash -c 'tox'

notebook-sanitizer: ## sanitize notebooks with configuration file
	@echo "Copying notebook output sanitizer ..."
	@-bash -c "curl -L $(SANITIZE_FILE) -o $(CURDIR)/docs/source/output-sanitize.cfg --silent"

test-notebooks: notebook-sanitizer  ## run notebook-based tests
	@echo "Running notebook-based tests"
	@$(MAKE) -f $(THIS_FILE) test-notebooks-impl

# Test one single notebook (add .run at the end of notebook path).
# Might require one time `make notebook-sanitizer`.
%.ipynb.run: %.ipynb
	@echo "Testing notebook $<"
	@$(MAKE) -f $(THIS_FILE) test-notebooks-impl NB_FILE="$<"

NB_FILE := $(CURDIR)/docs/source/notebooks/
test-notebooks-impl:
	@bash -c "env WPS_URL=$(WPS_URL) FINCH_WPS_URL=$(FINCH_WPS_URL) pytest --numprocesses=0 --nbval-lax --verbose $(NB_FILE) --sanitize-with $(CURDIR)/docs/source/output-sanitize.cfg --ignore $(CURDIR)/docs/source/notebooks/.ipynb_checkpoints"

ifeq "$(JUPYTER_NB_IP)" ""
JUPYTER_NB_IP := 0.0.0.0
endif
notebook: ## run Jupyter notebook server
	@echo "Running notebook server"
	@bash -c "env WPS_URL=$(WPS_URL) FINCH_WPS_URL=$(FINCH_WPS_URL) jupyter notebook --ip=$(JUPYTER_NB_IP) $(NB_FILE)"

coverage: ## check code coverage quickly with the default Python
	@bash -c 'coverage run --source raven -m pytest'
	@bash -c 'coverage report -m'
	@bash -c 'coverage html'
	$(BROWSER) htmlcov/index.html

# Only works for notebooks that passed ``make test-notebooks`` above.
# For those that failed, manually starting a local Jupyter server and refresh them manually.
refresh-notebooks: ## refreshing all notebook outputs under docs/source/notebooks
	@echo "Refresh all notebook outputs under docs/source/notebooks"
	@bash -c 'for nb in $(NB_FILE)/*.ipynb; do $(MAKE) -f $(THIS_FILE) refresh-notebooks-impl NB_REFRESH_FILE="$$nb"; done; cd $(APP_ROOT)'

# refresh one single notebook (add .refresh at the end of notebook path).
%.ipynb.refresh: %.ipynb
	@echo "Refreshing notebook $<"
	@$(MAKE) -f $(THIS_FILE) refresh-notebooks-impl NB_REFRESH_FILE="$<"

NB_REFRESH_FILE := ""
refresh-notebooks-impl: ## refresh one single notebook
	@bash -c 'WPS_URL="$(WPS_URL)" FINCH_WPS_URL="$(FINCH_WPS_URL)" jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=240 --output "$(NB_REFRESH_FILE)" --output-dir . "$(NB_REFRESH_FILE)"; sed -i "s@$(WPS_OUTPUT_URL)/@$(OUTPUT_URL)/@g" "$(NB_REFRESH_FILE)"'

## Sphinx targets:

docs: clean-docs ## generate Sphinx HTML documentation, including API docs
	@echo "Generating docs with Sphinx ..."
	@bash -c '$(MAKE) -C docs html'
	@echo "Open your browser to: file:/$(APP_ROOT)/docs/build/html/index.html"
	## do not execute xdg-open automatically since it hangs ReadTheDocs and job does not complete
	@echo "xdg-open $(APP_ROOT)/docs/build/html/index.html"

servedocs: docs ## compile the docs watching for changes
	@echo "Compiling the docs and watching for changes ..."
	@watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

## Docker targets:

docker-build: ## build the docker container
	@echo "Building the docker container ..."
	@docker build -t $(APP_NAME) .

docker-run: docker-build ## build and run the docker container locally
	@echo "Running the docker container locally ..."
	@docker run -d -p $(WPS_PORT):$(WPS_PORT) --name $(APP_NAME) $(APP_NAME)

## Deployment targets:

dist: clean ## build source and wheel package
	@python -m flit build
	@bash -c 'ls -l dist/'

release: dist ## upload source and wheel packages
	@echo "Uploading source and wheel packages ..."
	@python -m flit publish dist/*
