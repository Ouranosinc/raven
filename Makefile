# Application
APP_ROOT := $(CURDIR)
APP_NAME := raven

# Anaconda
ANACONDA_HOME ?= $(HOME)/miniconda
CONDA_ENV ?= $(APP_NAME)

# Choose Anaconda installer depending on your OS
ANACONDA_URL = https://repo.continuum.io/miniconda
RAVEN_URL = http://www.civil.uwaterloo.ca/jmai/raven/raven-rev163.zip
RAVEN_SRC = $(CURDIR)/src/RAVEN
UNAME_S := $(shell uname -s)
DOWNLOAD_CACHE = /tmp/

ifeq "$(UNAME_S)" "Linux"
FN := Miniconda3-latest-Linux-x86_64.sh
else ifeq "$(UNAME_S)" "Darwin"
FN := Miniconda3-latest-MacOSX-x86_64.sh
else
FN := unknown
endif

TEMP_FILES := *.egg-info *.log *.sqlite

# end of configuration

.DEFAULT_GOAL := all

.PHONY: all
all: help

.PHONY: help
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  help        to print this help message. (Default)"
	@echo "  install     to install $(APP_NAME) by running 'python setup.py develop'."
	@echo "  start       to start $(APP_NAME) service as daemon (background process)."
	@echo "  stop        to stop $(APP_NAME) service."
	@echo "  status      to show status of $(APP_NAME) service."
	@echo "  clean       to remove *all* files that are not controlled by 'git'. WARNING: use it *only* if you know what you do!"
	@echo "\nTesting targets:"
	@echo "  test        to run tests (but skip long running tests)."
	@echo "  testall     to run all tests (including long running tests)."
	@echo "  pep8        to run pep8 code style checks."
	@echo "\nSphinx targets:"
	@echo "  docs        to generate HTML documentation with Sphinx."

## Anaconda targets

.PHONY: anaconda
anaconda:
	@echo "Installing Anaconda ..."
	@test -d $(ANACONDA_HOME) || curl $(ANACONDA_URL)/$(FN) --silent --insecure --output "$(DOWNLOAD_CACHE)/$(FN)"
	@test -d $(ANACONDA_HOME) || bash "$(DOWNLOAD_CACHE)/$(FN)" -b -p $(ANACONDA_HOME)
	@echo "Please add '$(ANACONDA_HOME)/bin' to your PATH variable in '.bashrc'."

.PHONY: conda_env
conda_env: anaconda
	@echo "Updating conda environment $(CONDA_ENV) ..."
	"$(ANACONDA_HOME)/bin/conda" env update -n $(CONDA_ENV) -f environment.yml

## Build targets

.PHONY: bootstrap
bootstrap: conda_env bootstrap_dev
	@echo "Bootstrap ..."

.PHONY: bootstrap_dev
bootstrap_dev:
	@echo "Installing development requirements for tests and docs ..."
	@-bash -c "$(ANACONDA_HOME)/bin/conda install -y -n $(CONDA_ENV) pytest flake8 sphinx"
	@-bash -c "$(ANACONDA_HOME)/bin/conda install -c conda-forge -y -n $(CONDA_ENV) bumpversion"
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && pip install -r requirements_dev.txt"


.PHONY: raven_dev
raven_dev:
	@echo "Downloading RAVEN hydrological framework ..."
	@test -f $(CURDIR)/src/RAVEN.zip || curl $(RAVEN_URL) --output "$(CURDIR)/src/RAVEN.zip"
	@echo "Unzipping RAVEN ..."
	@test -d $(RAVEN_SRC) || unzip -j $(CURDIR)/src/RAVEN.zip -d "$(RAVEN_SRC)"
	@echo "Compiling RAVEN ..."
	@test -f $(RAVEN_SRC)/raven_rev.exe || $(MAKE) -C $(RAVEN_SRC) -j4
	@test -d bin || mkdir bin
	@-bash -c "cp $(RAVEN_SRC)/raven_rev.exe ./bin/raven"


.PHONY: raven_clean
raven_clean:
	@echo "Removing src and executable"
	@test -f $(CURDIR)/src/RAVEN.zip && rm -v "$(CURDIR)/src/RAVEN.zip" || echo "No zip to remove"
	@test -d $(RAVEN_SRC) && rm -rfv $(RAVEN_SRC) || echo "No src directory to remove"
	@test -f ./bin/raven && rm -v ./bin/raven || echo "No executable to remove"


.PHONY: install
install: bootstrap raven_dev
	@echo "Installing application ..."
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && python setup.py develop"
	@echo "\nStart service with \`make start'"

.PHONY: start
start:
	@echo "Starting application ..."
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && $(APP_NAME) start -d"

.PHONY: stop
stop:
	@echo "Stopping application ..."
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && $(APP_NAME) stop"

.PHONY: status
status:
	@echo "Show status ..."
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && $(APP_NAME) status"

.PHONY: clean
clean: srcclean envclean
	@echo "Cleaning generated files ..."
	@-for i in $(TEMP_FILES); do test -e $$i && rm -v -rf $$i; done

.PHONY: envclean
envclean:
	@echo "Removing conda env $(CONDA_ENV)"
	@-"$(ANACONDA_HOME)/bin/conda" remove -n $(CONDA_ENV) --yes --all

.PHONY: srcclean
srcclean:
	@echo "Removing *.pyc files ..."
	@-find $(APP_ROOT) -type f -name "*.pyc" -print | xargs rm

.PHONY: distclean
distclean: clean
	@echo "Cleaning ..."
	@git diff --quiet HEAD || echo "There are uncommited changes! Not doing 'git clean' ..."
	@-git clean -dfx -e *.bak -e custom.cfg

## Test targets

.PHONY: test
test:
	@echo "Running tests (skip slow and online tests) ..."
	@bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV);pytest -v -m 'not slow and not online'"

.PHONY: testall
testall:
	@echo "Running all tests (including slow and online tests) ..."
	@bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && pytest -v"

.PHONY: pep8
pep8:
	@echo "Running pep8 code style checks ..."
	@bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV) && flake8 raven tests"

##  Sphinx targets

.PHONY: docs
docs:
	@echo "Generating docs with Sphinx ..."
	@-bash -c "source $(ANACONDA_HOME)/bin/activate $(CONDA_ENV);$(MAKE) -C $@ clean html"
	@echo "open your browser: firefox docs/build/html/index.html"
