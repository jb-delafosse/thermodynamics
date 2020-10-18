import csv
import os

from thermodynamics.groupdecomposition.dto import PredefinedGroupDecompositionMethod

dirname = os.path.dirname(__file__)

PREDEFINED_DECOMPOSITION_TEST_FILE = {
    PredefinedGroupDecompositionMethod.UNIFAC_2013: os.path.join(
        dirname, "resources/UNIFAC_2013.csv"
    )
}


def pytest_generate_tests(metafunc):
    if "predefined_decomposition_method" in metafunc.fixturenames:
        for decomposition_method, path in PREDEFINED_DECOMPOSITION_TEST_FILE.items():
            with open(path, "r") as csvfile:
                csv_reader = csv.reader(csvfile, delimiter="\t")
                parameters = []
                for row in csv_reader:
                    parameters.append(
                        (decomposition_method, row[0], tuple((int(x) for x in row[1:])))
                    )
        metafunc.parametrize(
            "predefined_decomposition_method, smiles, expected_decomposition",
            parameters,
        )
