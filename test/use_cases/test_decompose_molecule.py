from typing import Generator

import pytest
from ordered_set import OrderedSet
from rdkit.Chem import AllChem, MolFromSmiles

from thermodynamics.groupdecomposition.dto import (
    FunctionalGroupDefinition,
    PredefinedGroupDecompositionMethod,
)
from thermodynamics.groupdecomposition.infrastructure.marshmallow_decomposition_loader import (
    MarshmallowPredefinedDecompositionLoader,
)
from thermodynamics.groupdecomposition.interfaces import (
    predefined_decomposition_loader as loader,
)
from thermodynamics.groupdecomposition.use_cases.decompose_molecule import (
    SimpleGroupDecomposition,
    load,
)


@pytest.fixture(autouse=True)
def mock_loader() -> Generator[loader.IPredefinedDecompositionLoader, None, None]:
    loader.PREDEFINED_DECOMPOSITION_LOADER = MarshmallowPredefinedDecompositionLoader()
    try:
        yield loader.PREDEFINED_DECOMPOSITION_LOADER
    finally:
        del loader.PREDEFINED_DECOMPOSITION_LOADER


def test_load():
    # GIVEN
    # - a molecule
    # - unifac decomposition
    molecule = AllChem.AddHs(MolFromSmiles("OCCCO"))
    unifac_decompo = load(PredefinedGroupDecompositionMethod.UNIFAC_2013)

    # WHEN
    # - using the decomposition on the molecule
    result = unifac_decompo.process((molecule,))

    # THEN
    # - the decomposition is performed correctly
    assert result


def test_process():
    # GIVEN
    # - Functional groups
    # - A molecule that can be decomposed entirely with the functional groups

    molecule = AllChem.AddHs(MolFromSmiles("OCCCO"))
    functional_groups = OrderedSet(
        (
            FunctionalGroupDefinition(smart="[OH1;$(*-C)][H]", name="OH", id=1),
            FunctionalGroupDefinition(smart="[H][CH2;!$(*=C)][H]", name="CH2", id=2),
        )
    )

    # WHEN
    # - Performing group decomposition on the molecule
    group_decompositions = SimpleGroupDecomposition(functional_groups).process(
        (molecule,)
    )

    # THEN
    # - The molecule is decomposed correctly
    assert group_decompositions[0][functional_groups[0]] == 2
    assert group_decompositions[0][functional_groups[1]] == 3


def test_process_not_entirely_decomposed():
    # GIVEN
    # - Functional groups
    # - A molecule that cannot be decomposed entirely with the functional groups
    molecule = AllChem.AddHs(MolFromSmiles("OCC(C)CO"))
    functional_groups = OrderedSet(
        (
            FunctionalGroupDefinition(smart="[OH1;$(*-C)][H]", name="OH", id=1),
            FunctionalGroupDefinition(smart="[H][CH2;!$(*=C)][H]", name="CH2", id=2),
        )
    )

    # WHEN
    # - Performing group decomposition on the molecule
    group_decompositions = SimpleGroupDecomposition(functional_groups).process(
        (molecule,)
    )

    # THEN
    # - The group_decomposition is None
    assert group_decompositions[0] is None
