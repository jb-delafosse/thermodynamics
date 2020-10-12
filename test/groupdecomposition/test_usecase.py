from unittest import TestCase

from ordered_set import OrderedSet
from rdkit.Chem import AllChem, MolFromSmiles

from thermodynamics.groupdecomposition.dto import FunctionalGroupDefinition
from thermodynamics.groupdecomposition.usecase import SimpleGroupDecomposition


class TestSimpleGroupDecomposition(TestCase):
    def test_process(self):
        # GIVEN
        # - Functional groups
        # - A molecule that can be decomposed entirely with the functional groups

        molecule = AllChem.AddHs(MolFromSmiles("OCCCO"))
        functional_groups = OrderedSet(
            (
                FunctionalGroupDefinition(smart="[OH1;$(*-C)][H]", name="OH"),
                FunctionalGroupDefinition(smart="[H][CH2;!$(*=C)][H]", name="CH2"),
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

    def test_process_not_entirely_decomposed(self):
        # GIVEN
        # - Functional groups
        # - A molecule that cannot be decomposed entirely with the functional groups
        molecule = AllChem.AddHs(MolFromSmiles("OCC(C)CO"))
        functional_groups = OrderedSet(
            (
                FunctionalGroupDefinition(smart="[OH1;$(*-C)][H]", name="OH"),
                FunctionalGroupDefinition(smart="[H][CH2;!$(*=C)][H]", name="CH2"),
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
