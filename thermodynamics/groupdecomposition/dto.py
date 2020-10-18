from dataclasses import dataclass
from enum import Enum
from typing import Tuple

from ordered_set import OrderedSet


@dataclass(frozen=True)
class FunctionalGroupDefinition:
    name: str
    smart: str
    id: int


GroupDecomposition = Tuple[int, ...]
DecompositionFunctionalGroups = OrderedSet[FunctionalGroupDefinition]


class PredefinedGroupDecompositionMethod(Enum):
    UNIFAC_2013 = 1
