from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict, tab

from .dispersion import DispersionModel
from .contact_partition import ContactPartitionModel
from .interfacial_area_density import InterfacialAreaDensityModel

@offbeat_define
class PairGeometryModel(OffbeatDict):
    dispersionModel: DispersionModel | None = None
    interfacialAreaDensityModel: InterfacialAreaDensityModel | None = None
    contactPartitionModel: ContactPartitionModel | None = None