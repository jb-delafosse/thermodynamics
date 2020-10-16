from unittest import TestCase

from thermodynamics.groupdecomposition.dto import PredefinedGroupDecompositionMethod
from thermodynamics.groupdecomposition.infrastructure.marshmallow_decomposition_loader import (
    MarshmallowPredefinedDecompositionLoader,
)


class TestMarshmallowPredefinedDecompositionLoader(TestCase):
    def test_load_UNIFAC_2013(self):
        # GIVEN / WHEN
        decomposition_functional_groups = MarshmallowPredefinedDecompositionLoader.load(
            PredefinedGroupDecompositionMethod.UNIFAC_2013
        )

        # THEN
        assert len(decomposition_functional_groups) == 153
