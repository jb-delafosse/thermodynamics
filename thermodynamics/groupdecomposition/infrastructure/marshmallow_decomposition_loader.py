import json
import os
from json import JSONDecodeError
from typing import Dict

from marshmallow import Schema, ValidationError, fields, post_load
from ordered_set import OrderedSet
from rdkit.Chem.rdmolfiles import MolFromSmarts

from thermodynamics.groupdecomposition.dto import (
    DecompositionFunctionalGroups,
    FunctionalGroupDefinition,
    PredefinedGroupDecompositionMethod,
)
from thermodynamics.groupdecomposition.interfaces.predefined_decomposition_loader import (
    IPredefinedDecompositionLoader,
)

dirname = os.path.dirname(__file__)

FGroupFilepath: Dict[PredefinedGroupDecompositionMethod, str] = {
    PredefinedGroupDecompositionMethod.UNIFAC_2013: os.path.join(
        dirname, "resources/UNIFAC_2013.json"
    )
}


class SmartStr(fields.Field):
    """Field that serializes to a SMART (string) and deserializes
    to a SMART (string). It performs validation on the SMART
    """

    def _validate(self, value):
        mol = MolFromSmarts(value)
        if not mol:
            raise ValidationError(f"Could not understand {value} as a SMART.")


class FunctionalGroupDefinitionSchema(Schema):
    smart = SmartStr()
    name = fields.Str()
    id = fields.Integer()

    @post_load
    def make_object(self, data, **kwargs):
        return FunctionalGroupDefinition(**data)


class MarshmallowPredefinedDecompositionLoader(IPredefinedDecompositionLoader):
    @staticmethod
    def load(
        decomposition_method: PredefinedGroupDecompositionMethod
    ) -> DecompositionFunctionalGroups:
        filepath = FGroupFilepath.get(decomposition_method)
        output = OrderedSet()
        if not filepath:
            raise NotImplementedError
        ids = []
        with open(filepath) as file:
            schema = FunctionalGroupDefinitionSchema()
            for line_count, line in enumerate(file):
                try:
                    # noinspection PyTypeChecker
                    data = json.loads(line)
                except JSONDecodeError as error:
                    raise ValueError(
                        f"Line {line_count} is not proper JSON."
                    ) from error
                try:
                    functional_group = schema.load(data)
                except (TypeError, ValidationError) as error:
                    raise ValueError(
                        f"Line {line_count} is not a proper functional group."
                    ) from error
                ids.append(functional_group.id)
                output.add(functional_group)
            if set(ids) != set(range(1, line_count + 2)):
                raise ValueError(
                    f"Error in Ids : Missing Index {set(range(1, line_count+2))-set(ids)}, Not expected index {set(ids)-set(range(1, line_count+2))}"
                )
        return output
