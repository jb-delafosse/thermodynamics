import json
import os
from json import JSONDecodeError

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

FGroupFilepath = {
    PredefinedGroupDecompositionMethod.UNIFAC_2013: os.path.join(
        dirname, "resources/UNIFAC_2013.json"
    )
}


class SmartStr(fields.Field):
    """Field that serializes to a SMART (string) and deserializes
    to a SMART (string). It performs validation on the SMART
    """

    def _serialize(self, value, attr, obj, **kwargs):
        return value

    def _deserialize(self, value, attr, data, **kwargs):
        mol = MolFromSmarts(value)
        if not mol:
            raise ValidationError(f"Could not understand {value} as a SMART.")
        return value


class FunctionalGroupDefinitionSchema(Schema):
    smart = SmartStr()
    name = fields.Str()

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
            raise NotImplementedError(
                f"{decomposition_method.name} is not implemented."
            )
        with open(filepath, "r") as file:
            schema = FunctionalGroupDefinitionSchema()
            for line_count, line in enumerate(file):
                try:
                    # noinspection PyTypeChecker
                    data = json.loads(line)
                except JSONDecodeError as error:
                    raise Exception(f"Line {line_count} is not proper JSON.") from error
                try:
                    functional_group = schema.load(data)
                except ValidationError as error:
                    raise Exception(
                        f"Line {line_count} is not a proper functional group."
                    ) from error

                output.add(functional_group)
        return output
