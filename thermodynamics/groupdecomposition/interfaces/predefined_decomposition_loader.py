from abc import ABC, abstractmethod

from thermodynamics.groupdecomposition.dto import (
    DecompositionFunctionalGroups,
    PredefinedGroupDecompositionMethod,
)


class IPredefinedDecompositionLoader(ABC):
    @staticmethod
    @abstractmethod
    def load(
        decomposition_method: PredefinedGroupDecompositionMethod
    ) -> DecompositionFunctionalGroups:
        """Load Functional Group Definition"""


PREDEFINED_DECOMPOSITION_LOADER: IPredefinedDecompositionLoader


def load_functional_groups(
    decomposition_method: PredefinedGroupDecompositionMethod
) -> DecompositionFunctionalGroups:
    return PREDEFINED_DECOMPOSITION_LOADER.load(decomposition_method)  # noqa: F821
