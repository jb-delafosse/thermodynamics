from dataclasses import dataclass
from typing import DefaultDict


@dataclass(frozen=True)
class FunctionalGroupDefinition:
    name: str
    smart: str


GroupDecompositionDict = DefaultDict[FunctionalGroupDefinition, int]
