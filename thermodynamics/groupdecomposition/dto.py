from dataclasses import dataclass
from enum import Enum
from typing import DefaultDict

from ordered_set import OrderedSet


@dataclass(frozen=True)
class FunctionalGroupDefinition:
    name: str
    smart: str


GroupDecompositionDict = DefaultDict[FunctionalGroupDefinition, int]
DecompositionFunctionalGroups = OrderedSet[FunctionalGroupDefinition]


class PredefinedGroupDecompositionMethod(Enum):
    UNIFAC_2013 = 1
