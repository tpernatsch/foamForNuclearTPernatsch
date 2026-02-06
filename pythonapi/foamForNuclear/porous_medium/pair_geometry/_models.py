from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict, tab

from .dispersion import DispersionModel
from .contact_partition import ContactPartitionModel
from .interfacial_area_density import InterfacialAreaDensityModel

@ffn_define
class PairGeometryModel(FoamForNuclearDict):
    dispersionModel: DispersionModel | None = None
    interfacialAreaDensityModel: InterfacialAreaDensityModel | None = None
    contactPartitionModel: ContactPartitionModel | None = None