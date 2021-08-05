import logging
from pathlib import Path
from typing import Any, List, Union

LOGGER = logging.getLogger("PYWPS")


def single_file_check(file_list: List[Union[str, Path]]) -> Any:
    """Return the first element of a file list. Raise an error if the list is empty or contains more than one element.

    Parameters
    ----------
    file_list : List[Union[str, Path]]
    """
    if isinstance(file_list, (str, Path)):
        return file_list

    try:
        if len(file_list) > 1:
            msg = "Multi-file handling for file is not supported. Exiting."
            raise NotImplementedError(msg)
        elif len(file_list) == 0:
            msg = "No files found. Exiting."
            raise FileNotFoundError(msg)
        return file_list[0]
    except (FileNotFoundError, NotImplementedError) as e:
        LOGGER.error(e)
        raise
