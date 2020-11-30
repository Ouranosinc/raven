"""Testing and tutorial utilities module."""
# Most of this code copied and adapted from xarray and xclim
import logging
from typing import Optional, Union
from pathlib import Path
from urllib.request import urlretrieve

from xarray import open_dataset as _open_dataset
from xarray import Dataset
from xarray.tutorial import file_md5_checksum

_default_cache_dir = Path.home() / ".raven_testing_data"

LOGGER = logging.getLogger("RAVEN")

__all__ = ["get_file", "open_dataset"]


def _get(
    fullname: Path, github_url: str, branch: str, suffix: str, cache_dir: Path,
) -> Path:
    cache_dir = cache_dir.absolute()
    local_file = cache_dir / fullname
    md5name = fullname.with_suffix("{}.md5".format(suffix))
    md5file = cache_dir / md5name

    if not local_file.is_file():
        # This will always leave this directory on disk.
        # We may want to add an option to remove it.
        local_file.parent.mkdir(parents=True, exist_ok=True)

        url = "/".join((github_url, "raw", branch, fullname.as_posix()))
        LOGGER.info("Fetching remote file: %s" % fullname.as_posix())
        urlretrieve(url, local_file)
        url = "/".join((github_url, "raw", branch, md5name.as_posix()))
        LOGGER.info("Fetching remote file md5: %s" % fullname.as_posix())
        urlretrieve(url, md5file)

        localmd5 = file_md5_checksum(local_file)
        try:
            with open(md5file) as f:
                remotemd5 = f.read()
            if localmd5 != remotemd5:
                local_file.unlink()
                msg = """
                    MD5 checksum does not match, try downloading dataset again.
                    """
                raise OSError(msg)
        except OSError as e:
            LOGGER.error(e)
    return local_file


# idea copied from xclim that borrowed it from xarray that was borrowed from Seaborn
def get_file(
    name,
    github_url: str = "https://github.com/Ouranosinc/raven-testdata",
    branch: Optional[str] = None,
    cache_dir: Path = _default_cache_dir,
) -> Path:
    """
    Open a dataset from the online GitHub-like repository.
    If a local copy is found then always use that to avoid network traffic.

    Parameters
    ----------
    name : str
        Name of the file containing the dataset including suffixes.
    github_url : str
        URL to Github repository where the data is stored.
    branch : str, optional
        For GitHub-hosted files, the branch to download from.
    cache_dir : Path
        The directory in which to search for and write cached data.

    Returns
    -------
    Union[Dataset, Path]

    See Also
    --------
    xarray.open_dataset
    """
    fullname = Path(name)
    suffix = fullname.suffix

    return _get(
        fullname=fullname,
        github_url=github_url,
        branch=branch,
        suffix=suffix,
        cache_dir=cache_dir,
    )


# idea copied from xclim that borrowed it from xarray that was borrowed from Seaborn
def open_dataset(
    name,
    suffix: Optional[str] = None,
    dap_url: Optional[str] = None,
    github_url: str = "https://github.com/Ouranosinc/raven-testdata",
    branch: Optional[str] = None,
    cache: bool = True,
    cache_dir: Path = _default_cache_dir,
    **kwds,
) -> Dataset:
    """
    Open a dataset from the online GitHub-like repository.
    If a local copy is found then always use that to avoid network traffic.

    Parameters
    ----------
    name : str
        Name of the file containing the dataset. If no suffix is given, assumed
        to be netCDF ('.nc' is appended).
    suffix : str, optional
        If no suffix is given, assumed to be netCDF ('.nc' is appended). For no suffix, set "".
    dap_url : str, optional
        URL to OPeNDAP folder where the data is stored. If supplied, supersedes github_url.
    github_url : str
        URL to Github repository where the data is stored.
    branch : str, optional
        For GitHub-hosted files, the branch to download from.
    cache_dir : Path
        The directory in which to search for and write cached data.
    cache : bool
        If True, then cache data locally for use on subsequent calls.
    kwds : dict, optional
        For NetCDF files, **kwds passed to xarray.open_dataset.

    Returns
    -------
    Union[Dataset, Path]

    See Also
    --------
    xarray.open_dataset
    """
    name = Path(name)
    if suffix is None:
        suffix = ".nc"
    fullname = name.with_suffix(suffix)

    if dap_url is not None:
        dap_file = Path(dap_url) / fullname
        try:
            ds = _open_dataset(dap_file, **kwds)
            return ds
        except OSError:
            msg = "OPeNDAP file not read. Verify that service is available."
            LOGGER.error(msg)
            raise

    local_file = _get(
        fullname=fullname,
        github_url=github_url,
        branch=branch,
        suffix=suffix,
        cache_dir=cache_dir,
    )

    try:
        ds = _open_dataset(local_file, **kwds)
        if not cache:
            ds = ds.load()
            local_file.unlink()
        return ds
    except OSError:
        raise
