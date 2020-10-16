from collections import defaultdict
from typing import DefaultDict, Optional, Tuple

from rdkit.Chem import Mol, MolFromSmarts

from thermodynamics.groupdecomposition.dto import (
    DecompositionFunctionalGroups,
    FunctionalGroupDefinition,
    GroupDecompositionDict,
    PredefinedGroupDecompositionMethod,
)
from thermodynamics.groupdecomposition.interfaces.predefined_decomposition_loader import (
    load_functional_groups,
)


class SimpleGroupDecomposition:
    def __init__(self, functional_groups: DecompositionFunctionalGroups):
        self._functional_groups = functional_groups

    def process(
        self, molecules: Tuple[Mol]
    ) -> Tuple[Optional[GroupDecompositionDict], ...]:
        return tuple(self._process_molecule(mol) for mol in molecules)

    def _process_molecule(
        self, molecule: Mol
    ) -> Optional[DefaultDict[FunctionalGroupDefinition, int]]:
        already_matched_atom = set()
        not_matched_atom = {atom.GetIdx() for atom in molecule.GetAtoms()}
        decomposition: GroupDecompositionDict = defaultdict(lambda: 0)
        for functional_group_def in self._functional_groups:
            if not not_matched_atom:
                break
            for match in molecule.GetSubstructMatches(
                MolFromSmarts(functional_group_def.smart)
            ):
                if all(atom in not_matched_atom for atom in match):
                    not_matched_atom = not_matched_atom - set(match)
                    already_matched_atom.update(set(match))
                    decomposition[functional_group_def] += 1
        if not_matched_atom:
            return None
        else:
            return decomposition


def load(
    predefined_decomposition: PredefinedGroupDecompositionMethod
) -> SimpleGroupDecomposition:
    functional_groups = load_functional_groups(predefined_decomposition)
    return SimpleGroupDecomposition(functional_groups)
