from unittest.mock import mock_open, patch

import pytest

from thermodynamics.groupdecomposition.dto import PredefinedGroupDecompositionMethod
from thermodynamics.groupdecomposition.infrastructure.marshmallow_decomposition_loader import (
    MarshmallowPredefinedDecompositionLoader,
)


@pytest.mark.parametrize(
    "line, expected_exception",
    [
        ("not json", ValueError),
        ('{"name": "OH"}', ValueError),
        ('{"name": "OH", "smart": "not-a-smart"}', ValueError),
    ],
)
def test_eval(line, expected_exception):
    # GIVEN
    # - a faulty line in the jsonlines file
    # WHEN
    # - trying to load the faulty file
    # THEN
    # - the expected exception is raised
    with patch("builtins.open", new_callable=mock_open, read_data=line):
        with pytest.raises(expected_exception):
            MarshmallowPredefinedDecompositionLoader.load(
                PredefinedGroupDecompositionMethod.UNIFAC_2013
            )


def test_load_unifac_2013():
    # GIVEN
    # - the unifac 2013 decomposition file
    # WHEN
    # - loading the function groups
    decomposition_functional_groups = MarshmallowPredefinedDecompositionLoader.load(
        PredefinedGroupDecompositionMethod.UNIFAC_2013
    )

    # THEN
    # - the 153 groups load correctly
    assert len(decomposition_functional_groups) == 153
