from typing import List, Optional

import numpy as np


def get_r_and_q(
    subgroups: Optional[List[str]], nu: np.ndarray, smiles: Optional[str] = ""
):
    """
    Calculates van der Waals volume (r_i) and surface area (q_i) for a single
        component

    Incoming parameters:
        subgroups = string list, subgroups present in molecule
        nu = array, number of occurences for each subgroup,
            in the same order as the "subgroups" list
    Optional parameter:
        smiles = string, determine the subgroup and occurrence automatically
            from molecule SMILES (***not available yet***)
    Outgoing parameters:
        r = van der Waals volume for component i
        q = van der Waals surface area for component i
        R_pure = van der Waals volume for each subgroup in component i
        Q_pure = van der Waals surface area for each subgroup in component i
    Prerequisite:
        File "unifac_R_Q_params.txt" must be present in same folder as
            module
    Algorithm for subgroups and nu:
        - If subgroups are provided (most certainly nu will be too), then
            the parameter "smiles" is ignored, and use provided subgroup
            and nus by default
        - To request automatic decomposition using SMILES, do the following:
            1. Send empty lists for subgroup and nu (i.e., [], np.array([]))
            2. Provide SMILES as a string
        - If SMILES is not provided, then subgroups and nu must be provided,
            and vice-versa
        - If all of subgroups, nu and smile are provided, the function uses
            the provided subgroup and nu from the user
    Notes regarding "unifac_R_Q_params.txt":
        1. The first 6 lines are explanations and have no use in the code
        2. The 35th word (excl. spaces) marks the beginning of the parameters
        3. Subgroup names are unique
        4. The entry following the subgroup name will be its R value
        5. The second entry following the subgroup name will be its Q value
    Reference:
        J. Gmehling, B. Kolbe, M. Kleiber, J. Rarey, Chemical Thermodynamics
            for Process Simulation, Wiley-VCH, Second Edition, 2019
    """

    filename = "unifac_R_Q_params.txt"  # Modify file location if needed

    if not subgroups:  # No subgroups (and nu) provided
        raise NotImplementedError(
            "Automatic Decomposition not implemented."
        )  # TODO: call automatic decomposition

    if not subgroups and not smiles:
        raise Exception(
            "If Subgroups are missing, must provide SMILES for automatic decomposition."
        )

    with open(filename) as f_obj:
        params_str = f_obj.read()  # Read all contents

    params = params_str.split()  # Remove all spacing
    params = params[35:]  # Remove all unnecessary comments

    R_pure = np.empty(len(subgroups))
    Q_pure = np.empty(len(subgroups))

    for k in range(len(subgroups)):
        index = params.index(subgroups[k])  # Locate subgroup line in file
        R_pure[k] = params[index + 1]
        Q_pure[k] = params[index + 2]

    r = sum(nu * R_pure)  # Eq.(5.87)
    q = sum(nu * Q_pure)  # Eq.(5.88)

    return r, q, R_pure, Q_pure


def combinatorial(r, q, x, Z=10):
    """
    Calculates the combinatorial part of the activity coefficient:
    ln(gamma_i^C)

    Mandatory parameter:
        r = float array, van der Waals volumes for n components
        q = float array, van der Waals surface areas for n components
        x = float array, mole fraction for n components
    Optional parameter:
        Z = integer, coordination number of the combinatorial term. Default
            value is 10
    Return:
        lngammaC = float array, natural log of the combinatorial part of the
            activity coefficient for n components
    Reference:
        J. Gmehling, B. Kolbe, M. Kleiber, J. Rarey, Chemical Thermodynamics
            for Process Simulation, Wiley-VCH, Second Edition, 2019
    """

    V = r / sum(r * x)  # Vi (volume/mole fraction ratio), Eq. (5.85)

    F = q / sum(q * x)  # Fi (surface area/mole fraction ratio), Eq. (5.86)

    lngammaC = 1 - V + np.log(V) - 5 * q * (1 - (V / F) + np.log(V / F))

    return lngammaC
