from pathlib import Path

import pytest

from raven.utilities import single_file_check


class TestUtilities:
    def test_single_file_check(self):
        one = [Path(__file__).parent / "__init__.py"]
        zero = list()
        three = [1, Path().root, 2.333]

        assert single_file_check(one) == one[0]

        with pytest.raises(FileNotFoundError):
            single_file_check(zero)

        with pytest.raises(NotImplementedError):
            single_file_check(three)
