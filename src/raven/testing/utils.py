"""Tools for searching for and acquiring test data."""

from __future__ import annotations

import importlib.metadata as ilm
import importlib.resources as ilr
import logging
import os
import re
import tempfile
import time
import warnings
from collections.abc import Callable
from datetime import datetime as dt
from functools import wraps
from io import StringIO
from pathlib import Path
from shutil import copytree
from typing import IO, Any, TextIO
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import urlretrieve

from filelock import FileLock
from packaging.version import Version
from xarray import Dataset
from xarray import open_dataset as _open_dataset
from xclim.testing.utils import show_versions as _show_versions

import raven

try:
    import pooch
except ImportError:
    warnings.warn(
        "The `pooch` library is not installed. The default cache directory for testing data will not be set."
    )
    pooch = None

LOGGER = logging.getLogger("raven.testing.utils")

__all__ = [
    "TESTDATA_BRANCH",
    "TESTDATA_CACHE_DIR",
    "TESTDATA_REPO_URL",
    "audit_url",
    "default_testdata_cache",
    "default_testdata_repo_url",
    "default_testdata_version",
    "gather_testing_data",
    "open_dataset",
    "populate_testing_data",
    "show_versions",
    "testing_setup_warnings",
    "yangtze",
]

default_testdata_version = "v2025.6.12"
"""Default version of the testing data to use when fetching datasets."""

default_testdata_repo_url = (
    "https://raw.githubusercontent.com/Ouranosinc/raven-testdata/"
)
"""Default URL of the testing data repository to use when fetching datasets."""

try:
    default_testdata_cache = Path(pooch.os_cache("raven-testdata"))
    """Default location for the testing data cache."""
except (AttributeError, TypeError):
    default_testdata_cache = None

TESTDATA_REPO_URL = str(os.getenv("RAVEN_TESTDATA_REPO_URL", default_testdata_repo_url))
"""
Sets the URL of the testing data repository to use when fetching datasets.

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export RAVEN_TESTDATA_REPO_URL="https://github.com/my_username/raven-testdata"

or setting the variable at runtime:

.. code-block:: console

    $ env RAVEN_TESTDATA_REPO_URL="https://github.com/my_username/raven-testdata" pytest
"""

TESTDATA_BRANCH = str(os.getenv("RAVEN_TESTDATA_BRANCH", default_testdata_version))
"""
Sets the branch of the testing data repository to use when fetching datasets.

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export RAVEN_TESTDATA_BRANCH="my_testing_branch"

or setting the variable at runtime:

.. code-block:: console

    $ env RAVEN_TESTDATA_BRANCH="my_testing_branch" pytest
"""

TESTDATA_CACHE_DIR = os.getenv("RAVEN_TESTDATA_CACHE_DIR", default_testdata_cache)
"""
Sets the directory to store the testing datasets.

If not set, the default location will be used (based on ``platformdirs``, see :func:`pooch.os_cache`).

Notes
-----
When running tests locally, this can be set for both `pytest` and `tox` by exporting the variable:

.. code-block:: console

    $ export RAVEN_TESTDATA_CACHE_DIR="/path/to/my/data"

or setting the variable at runtime:

.. code-block:: console

    $ env RAVEN_TESTDATA_CACHE_DIR="/path/to/my/data" pytest
"""


def show_versions(
    file: os.PathLike | StringIO | TextIO | None = None,
    deps: list | None = None,
) -> str | None:
    """
    Print the versions of RavenWPS and its dependencies.

    Parameters
    ----------
    file : {os.PathLike, StringIO, TextIO}, optional
        If provided, prints to the given file-like object. Otherwise, returns a string.
    deps : list, optional
        A list of dependencies to gather and print version information from. Otherwise, print RavenWPS dependencies.

    Returns
    -------
    str or None
        The formatted version information if `file` is not provided, otherwise None.
    """

    def _get_raven_dependencies():
        raven_metadata = ilm.metadata("raven")
        requires = raven_metadata.get_all("Requires-Dist")
        requires = [
            re.match(r"^[A-Za-z0-9_.\-]+", req).group(0)
            for req in requires
            if re.match(r"^[A-Za-z0-9_.\-]+", req)
        ]
        sorted_deps = sorted(list(set(requires) - {"raven"}))

        return ["raven"] + sorted_deps

    if deps is None:
        deps = _get_raven_dependencies()

    return _show_versions(file=file, deps=deps)


# Test Data Utilities ###


def testing_setup_warnings():
    """Warn users about potential incompatibilities between RavenWPS and raven-testdata versions."""
    if (
        re.match(r"^\d+\.\d+\.\d+$", raven.__version__)
        and TESTDATA_BRANCH != default_testdata_version
    ):
        # This does not need to be emitted on GitHub Workflows and ReadTheDocs
        if not os.getenv("CI") and not os.getenv("READTHEDOCS"):
            warnings.warn(
                f"`birdhouse-raven` stable ({raven.__version__}) is running tests against a non-default "
                f"branch of the testing data. It is possible that changes to the testing data may "
                f"be incompatible with some assertions in this version. "
                f"Please be sure to check {TESTDATA_REPO_URL} for more information.",
            )

    if re.match(r"^v\d+\.\d+\.\d+", TESTDATA_BRANCH):
        # Find the date of the last modification of RavenWPS source files to generate a calendar version
        install_date = dt.strptime(
            time.ctime(Path(raven.__file__).stat().st_mtime),
            "%a %b %d %H:%M:%S %Y",
        )
        install_calendar_version = (
            f"{install_date.year}.{install_date.month}.{install_date.day}"
        )

        if Version(TESTDATA_BRANCH) > Version(install_calendar_version):
            warnings.warn(
                f"The installation date of `birdhouse-raven` ({install_date.ctime()}) "
                f"predates the last release of testing data ({TESTDATA_BRANCH}). "
                "It is very likely that the testing data is incompatible with this build of RavenWPS.",
            )


def load_registry(
    branch: str = TESTDATA_BRANCH,
    repo: str = TESTDATA_REPO_URL,
    force_download: bool = False,
) -> dict[str, str]:
    """
    Load the registry file for the test data.

    Parameters
    ----------
    branch : str
        Branch of the repository to use when fetching testing datasets.
    repo : str
        URL of the repository to use when fetching testing datasets.
    force_download : bool
        If True, force the download of the registry file even if it already exists.

    Returns
    -------
    dict
        Dictionary of filenames and hashes.
    """

    def load_registry_from_file(
        _registry_file: str | Path,
    ) -> dict[str, str]:
        """Load the registry from a file."""
        with Path(_registry_file).open(encoding="utf-8") as f:
            return {line.split()[0]: line.split()[1] for line in f}

    if not repo.endswith("/"):
        repo = f"{repo}/"
    remote_registry = audit_url(
        urljoin(
            urljoin(repo, branch if branch.endswith("/") else f"{branch}/"),
            "data/registry.txt",
        )
    )

    if repo != default_testdata_repo_url:
        external_repo_name = urlparse(repo).path.split("/")[-2]
        external_branch_name = branch.split("/")[-1]
        testing_folder = Path(str(ilr.files("raven").joinpath("testing")))
        registry_file = testing_folder.joinpath(
            f"registry.{external_repo_name}.{external_branch_name}.txt"
        )
        lockfile = testing_folder.joinpath(".lock")
        with FileLock(lockfile):
            if not registry_file.exists():
                urlretrieve(remote_registry, registry_file)  # noqa: S310
        lockfile.unlink(missing_ok=True)

    elif branch != default_testdata_version:
        # If force_download is True, download to a transient directory for testing purposes
        if force_download:
            with tempfile.TemporaryDirectory() as tmp_dir:
                custom_registry_folder = Path(tmp_dir).joinpath("testing", branch)
                custom_registry_folder.mkdir(parents=True, exist_ok=True)
                registry_file = custom_registry_folder.joinpath("registry.txt")
                urlretrieve(remote_registry, registry_file)  # noqa: S310
                return load_registry_from_file(registry_file)
        else:
            # If the branch is not the default version, check if the registry file exists
            custom_registry_folder = Path(
                str(ilr.files("raven").joinpath("testing", branch))
            )
            custom_registry_folder.mkdir(parents=True, exist_ok=True)
            registry_file = custom_registry_folder.joinpath("registry.txt")
            with FileLock(custom_registry_folder.joinpath(".lock")):
                if not registry_file.exists():
                    urlretrieve(remote_registry, registry_file)  # noqa: S310
            return load_registry_from_file(registry_file)

    else:
        registry_file = Path(str(ilr.files("raven").joinpath("testing/registry.txt")))

    if not registry_file.exists():
        msg = f"Registry file not found: {registry_file}"
        raise FileNotFoundError(msg)

    return load_registry_from_file(registry_file)


def yangtze(
    repo: str = TESTDATA_REPO_URL,
    branch: str = TESTDATA_BRANCH,
    cache_dir: str | Path = TESTDATA_CACHE_DIR,
    allow_updates: bool = True,
    force_download: bool = False,
):
    """
    Pooch registry instance for RavenWPS test data.

    Parameters
    ----------
    repo : str
        URL of the repository to use when fetching testing datasets.
    branch : str
        Branch of repository to use when fetching testing datasets.
    cache_dir : str or Path
        The path to the directory where the data files are stored.
    allow_updates : bool
        If True, allow updates to the data files. Default is True.
    force_download : bool
        If True, force the download of the registry file even if it already exists.

    Returns
    -------
    pooch.Pooch
        The Pooch instance for accessing the RavenWPS testing data.

    Notes
    -----
    There are three environment variables that can be used to control the behaviour of this registry:
        - ``RAVEN_TESTDATA_CACHE_DIR``: If this environment variable is set, it will be used as the
          base directory to store the data files.
          The directory should be an absolute path (i.e. it should start with ``/``).
          Otherwise, the default location will be used (based on ``platformdirs``, see :py:func:`pooch.os_cache`).
        - ``RAVEN_TESTDATA_REPO_URL``: If this environment variable is set, it will be used as the URL of
          the repository to use when fetching datasets. Otherwise, the default repository will be used.
        - ``RAVEN_TESTDATA_BRANCH``: If this environment variable is set, it will be used as the branch of
          the repository to use when fetching datasets. Otherwise, the default branch will be used.

    Examples
    --------
    Using the registry to download a file:

    .. code-block:: python

        import xarray as xr
        from raven.testing import yangtze

        example_file = yangtze().fetch("example.nc")
        data = xr.open_dataset(example_file)
    """
    if pooch is None:
        raise ImportError(
            "The `pooch` package is required to fetch the RavenWPS testing data. "
            "You can install it with `pip install pooch` or `pip install birdhouse-raven[dev]`."
        )
    if not repo.endswith("/"):
        repo = f"{repo}/"
    remote = audit_url(
        urljoin(urljoin(repo, branch if branch.endswith("/") else f"{branch}/"), "data")
    )

    _yangtze = pooch.create(
        path=cache_dir,
        base_url=remote,
        version=default_testdata_version,
        version_dev=branch,
        allow_updates=allow_updates,
        registry=load_registry(branch=branch, repo=repo, force_download=force_download),
    )

    # Add a custom fetch method to the Pooch instance
    # Needed to address: https://github.com/readthedocs/readthedocs.org/issues/11763
    # Fix inspired by @bjlittle (https://github.com/bjlittle/geovista/pull/1202)
    _yangtze.fetch_diversion = _yangtze.fetch

    # Overload the fetch method to add user-agent headers
    @wraps(_yangtze.fetch_diversion)
    def _fetch(
        *args, **kwargs: bool | Callable
    ) -> str:  # numpydoc ignore=GL08  # *args: str
        def _downloader(
            url: str,
            output_file: str | IO,
            poocher: pooch.Pooch,
            check_only: bool | None = False,
        ) -> None:
            """Download the file from the URL and save it to the save_path."""
            headers = {"User-Agent": f"RavenWPS ({raven.__version__})"}
            downloader = pooch.HTTPDownloader(headers=headers)
            return downloader(url, output_file, poocher, check_only=check_only)

        # default to our http/s downloader with user-agent headers
        kwargs.setdefault("downloader", _downloader)
        return _yangtze.fetch_diversion(*args, **kwargs)

    # Replace the fetch method with the custom fetch method
    _yangtze.fetch = _fetch

    return _yangtze


def open_dataset(
    name: str,
    _yangtze_kwargs: dict[str, Path | str | bool] | None = None,
    **xr_kwargs: Any,
) -> Dataset:
    r"""
    Convenience function to open a dataset from the RavenWPS testing data using the `yangtze` class.

    This is a thin wrapper around the `yangtze` class to make it easier to open RavenWPS testing datasets.

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
    _yangtze_kwargs : dict
        Keyword arguments passed to the yangtze function.
    **xr_kwargs : Any
        Keyword arguments passed to xarray.open_dataset.

    Returns
    -------
    xarray.Dataset
        The dataset.

    See Also
    --------
    xarray.open_dataset : Open and read a dataset from a file or file-like object.
    yangtze : Pooch wrapper for accessing the RavenWPS testing data.
    """
    if _yangtze_kwargs is None:
        _yangtze_kwargs = {}
    return _open_dataset(yangtze(**_yangtze_kwargs).fetch(name), **xr_kwargs)


def populate_testing_data(
    temp_folder: Path | None = None,
    repo: str = TESTDATA_REPO_URL,
    branch: str = TESTDATA_BRANCH,
    retry: int = 3,
    local_cache: Path = TESTDATA_CACHE_DIR,
) -> None:
    """
    Populate the local cache with the testing data.

    Parameters
    ----------
    temp_folder : Path, optional
        Path to a temporary folder to use as the local cache. If not provided, the default location will be used.
    repo : str, optional
        URL of the repository to use when fetching testing datasets.
    branch : str, optional
        Branch of raven-testdata to use when fetching testing datasets.
    retry : int
        Number of times to retry downloading the files in case of failure. Default: 3.
    local_cache : Path
        The path to the local cache. Defaults to the location set by the platformdirs library.
        The testing data will be downloaded to this local cache.
    """
    # Create the Pooch instance
    n = yangtze(repo=repo, branch=branch, cache_dir=temp_folder or local_cache)

    # Download the files
    errored_files = []
    for file in load_registry():
        msg = f"Downloading file `{file}` from remote repository..."
        logging.info(msg)
        for attempt in range(retry):
            try:
                n.fetch(file)
            except HTTPError:  # noqa: PERF203
                msg = f"Failed to download file `{file}` on attempt {attempt + 1}."
                logging.info(msg)
            else:
                logging.info("File was downloaded successfully.")
                break
        else:
            msg = f"Failed to download file `{file}` after {retry} attempts."
            logging.error(msg)
            errored_files.append(file)

    if errored_files:
        msg = f"The following files were unable to be downloaded: {errored_files}"
        logging.error(msg)


def gather_testing_data(
    worker_cache_dir: str | os.PathLike[str] | Path,
    worker_id: str,
    _cache_dir: str | os.PathLike[str] | None = TESTDATA_CACHE_DIR,
) -> None:
    """
    Gather testing data across workers.

    Parameters
    ----------
    worker_cache_dir : str or Path
        The directory to store the testing data.
    worker_id : str
        The worker ID.
    _cache_dir : str or Path, optional
        The directory to store the testing data. Default is None.

    Raises
    ------
    ValueError
        If the cache directory is not set.
    FileNotFoundError
        If the testing data is not found.
    """
    if _cache_dir is None:
        raise ValueError(
            "The cache directory must be set. Please set the `cache_dir` parameter or the `RAVEN_TESTDATA_CACHE_DIR` environment variable."
        )
    cache_dir = Path(_cache_dir)

    if worker_id == "master":
        populate_testing_data(branch=TESTDATA_BRANCH)
    else:
        cache_dir.mkdir(exist_ok=True, parents=True)
        lockfile = cache_dir.joinpath(".lock")
        test_data_being_written = FileLock(lockfile)
        with test_data_being_written:
            # This flag prevents multiple calls from re-attempting to download testing data in the same pytest run
            populate_testing_data(branch=TESTDATA_BRANCH)
            cache_dir.joinpath(".data_written").touch()
        with test_data_being_written.acquire():
            if lockfile.exists():
                lockfile.unlink()
        copytree(cache_dir.joinpath(default_testdata_version), worker_cache_dir)


# Testing Utilities ###


def audit_url(url: str, context: str | None = None) -> str:
    """
    Check if the URL is well-formed.

    Parameters
    ----------
    url : str
        The URL to check.
    context : str, optional
        Additional context to include in the error message. Default is None.

    Returns
    -------
    str
        The URL if it is well-formed.

    Raises
    ------
    URLError
        If the URL is not well-formed.
    """
    msg = ""
    result = urlparse(url)
    if result.scheme == "http":
        msg = f"{context if context else ''} URL is not using secure HTTP: '{url}'".strip()
    if not all([result.scheme, result.netloc]):
        msg = f"{context if context else ''} URL is not well-formed: '{url}'".strip()

    if msg:
        LOGGER.error(msg)
        raise URLError(msg)
    return url
