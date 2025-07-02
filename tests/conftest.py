from __future__ import annotations

import logging
from pathlib import Path

import pytest

from raven.testing.utils import (
    TESTDATA_BRANCH,
    TESTDATA_CACHE_DIR,
    TESTDATA_REPO_URL,
    default_testdata_cache,
    gather_testing_data,
)
from raven.testing.utils import (
    testing_setup_warnings,
)
from raven.testing.utils import yangtze as _yangtze


@pytest.fixture(scope="session")
def threadsafe_data_dir(tmp_path_factory) -> Path:
    return Path(tmp_path_factory.getbasetemp().joinpath("data"))


@pytest.fixture(scope="session")
def yangtze(threadsafe_data_dir, worker_id):
    return _yangtze(
        repo=TESTDATA_REPO_URL,
        branch=TESTDATA_BRANCH,
        cache_dir=(
            TESTDATA_CACHE_DIR if worker_id == "master" else threadsafe_data_dir
        ),
    )


@pytest.fixture(autouse=True, scope="session")
def gather_session_data(request, yangtze, worker_id):
    """
    Gather testing data on pytest run.

    When running pytest with multiple workers, one worker will copy data remotely to the default cache dir while
    other workers wait using a lockfile. Once the lock is released, all workers will then copy data to their local
    threadsafe_data_dir. As this fixture is scoped to the session, it will only run once per pytest run.
    """
    testing_setup_warnings()
    gather_testing_data(worker_cache_dir=yangtze.path, worker_id=worker_id)

    def remove_data_written_flag():
        """Clean up the cache folder once we are finished."""
        flag = default_testdata_cache.joinpath(".data_written")
        if flag.exists():
            try:
                flag.unlink()
            except FileNotFoundError:
                logging.info(
                    "Teardown race condition occurred: .data_written flag already removed. Lucky!"
                )
                pass

    request.addfinalizer(remove_data_written_flag)
