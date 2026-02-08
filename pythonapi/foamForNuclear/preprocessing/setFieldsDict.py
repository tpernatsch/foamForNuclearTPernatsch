from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import *
from foamForNuclear.fields import Field
from foamForNuclear.openfoamFile import OpenFOAMFile


def getFieldValuesAsText(fieldValues: list[tuple[Field, float]], depth: int=0):
    text = f"\n{depth*tab}(\n"
    for field, value in fieldValues:
        text += f"{(depth+1)*tab}{field.fieldType}Value {field.name} {value}\n"
    text += f"{depth*tab})"
    return(text)


class SetFieldRegion(OpenFOAMDict):
    """
    Base class for set values in cells.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.

    Attributes
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]]
        ):
        super().__init__()

        for field, value in fieldValues:
            check_type("field", field, Field)
            check_type("value", value, (float, int, Vector))

        self.fieldValues: list[tuple[Field, float]] = fieldValues

    def __repr__(self, depth: int=0):
        self.__setitem__("fieldValues", getFieldValuesAsText(self.fieldValues, depth=depth+1))
        return super().__repr__(depth)


    def add_field_values(self, field: Field, value: float | Vector):
        check_type("field", field, Field)
        check_type("value", value, (int, float, Vector))
        self.fieldValues.append((field, value))


class CylinderToCell(SetFieldRegion):
    """
    Set values in cells based on a cylinder.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    radius : float
        Radius of the cylinder.
    point1 : Vector
        Center of the first face of the cylinder.
    point2 : Vector
        Center of the second face of the cylinder.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]],
            radius: float,
            point1: Vector,
            point2: Vector,
        ):
        super().__init__(fieldValues=fieldValues)

        self.name = "cylinderToCell"

        self.radius = radius
        self.point1 = point1
        self.point2 = point2

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius) -> None:
        check_type("radius", radius, (int, float))
        self._radius = radius
        self.__setitem__('radius', radius)

    @property
    def point1(self):
        return self._point1

    @point1.setter
    def point1(self, point1) -> None:
        check_type("point1", point1, Vector)
        self._point1 = point1
        self.__setitem__('point1', point1)

    @property
    def point2(self):
        return self._point2

    @point2.setter
    def point2(self, point2) -> None:
        check_type("point2", point2, Vector)
        self._point2 = point2
        self.__setitem__('point2', point2)


class CylinderAnnulusToCell(SetFieldRegion):
    """
    Set values in cells based on a hollow cylinder.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    innerRadius : float
        Inner radius of the cylinder.
    outerRadius : float
        Outer radius of the cylinder.
    point1 : Vector
        Center of the first face of the cylinder.
    point2 : Vector
        Center of the second face of the cylinder.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]],
            innerRadius: float,
            outerRadius: float,
            point1: Vector,
            point2: Vector,
        ):
        super().__init__(fieldValues=fieldValues)

        self.name = "cylinderAnnulusToCell"

        self.innerRadius = innerRadius
        self.outerRadius = outerRadius
        self.point1 = point1
        self.point2 = point2

    @property
    def outerRadius(self):
        return self._outerRadius

    @outerRadius.setter
    def outerRadius(self, outerRadius) -> None:
        check_type("outerRadius", outerRadius, (int, float))
        self._outerRadius = outerRadius
        self.__setitem__('radius', outerRadius)

    @property
    def innerRadius(self):
        return self._innerRadius

    @innerRadius.setter
    def innerRadius(self, innerRadius) -> None:
        check_type("innerRadius", innerRadius, (int, float))
        self._innerRadius = innerRadius
        self.__setitem__('innerRadius', innerRadius)

    @property
    def point1(self):
        return self._point1

    @point1.setter
    def point1(self, point1) -> None:
        check_type("point1", point1, Vector)
        self._point1 = point1
        self.__setitem__('point1', point1)

    @property
    def point2(self):
        return self._point2

    @point2.setter
    def point2(self, point2) -> None:
        check_type("point2", point2, Vector)
        self._point2 = point2
        self.__setitem__('point2', point2)


class SphereToCell(SetFieldRegion):
    """
    Set values in cells based on a hollow sphere.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    innerRadius : float
        Inner radius of the sphere.
    outerRadius : float
        Outer radius of the sphere.
    origin : Vector
        Center of the sphere.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]],
            innerRadius: float,
            outerRadius: float,
            origin: Vector,
        ):
        super().__init__(fieldValues=fieldValues)

        self.name = "sphereToCell"

        self.innerRadius = innerRadius
        self.outerRadius = outerRadius
        self.origin = origin

    @property
    def outerRadius(self):
        return self._outerRadius

    @outerRadius.setter
    def outerRadius(self, outerRadius) -> None:
        check_type("outerRadius", outerRadius, (int, float))
        self._outerRadius = outerRadius
        self.__setitem__('radius', outerRadius)

    @property
    def innerRadius(self):
        return self._innerRadius

    @innerRadius.setter
    def innerRadius(self, innerRadius) -> None:
        check_type("innerRadius", innerRadius, (int, float))
        self._innerRadius = innerRadius
        self.__setitem__('innerRadius', innerRadius)

    @property
    def point1(self):
        return self._point1

    @point1.setter
    def point1(self, point1) -> None:
        check_type("point1", point1, Vector)
        self._point1 = point1
        self.__setitem__('point1', point1)

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin) -> None:
        check_type("origin", origin, Vector)
        self._origin = origin
        self.__setitem__('origin', origin)


class BoxToCell(SetFieldRegion):
    """
    Set values in cells based on a box.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    lowCorner : Vector
        Coordinates of the lower corner of the parallelepiped.
    highCorner : Vector
        Coordinates of the upper corner of the parallelepiped.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]],
            lowCorner: Vector,
            highCorner: Vector,
        ):
        super().__init__(fieldValues=fieldValues)

        self.name = "boxToCell"

        self.box: Box = Box(lowCorner, highCorner)

    def __repr__(self, depth=0):
        self.__setitem__("box", self.box)
        return super().__repr__(depth)

    @property
    def lowCorner(self):
        return self._lowCorner

    @lowCorner.setter
    def lowCorner(self, lowCorner) -> None:
        check_type("lowCorner", lowCorner, Vector)
        self._box.lowCornder = lowCorner

    @property
    def highCorner(self):
        return self._highCorner

    @highCorner.setter
    def highCorner(self, highCorner) -> None:
        check_type("highCorner", highCorner, Vector)
        self._box.lowCornder = highCorner


class ZoneToCell(SetFieldRegion):
    """
    Set values in cells based on a cellZone.

    Parameters
    ----------
    fieldValues : list[tuple[Field, float]]
        List of field with the corresponding value to assign.
    zoneName : str
        Name of the cellZone.
    """
    def __init__(
            self,
            fieldValues: list[tuple[Field, float]],
            zoneName: str,
        ):
        super().__init__(fieldValues)

        self.name = "zoneToCell"

        self.zoneName = zoneName

    @property
    def zoneName(self):
        return self._zoneName

    @zoneName.setter
    def zoneName(self, zoneName) -> None:
        check_type("zoneName", zoneName, str)
        self._zoneName = zoneName
        self.__setitem__('zone', zoneName)


class SetFieldsDict(OpenFOAMFile):
    """
    Set values on a selected set of cells/patch-faces via a dictionary.

    Parameters
    ----------
    region : str
        Name of the region.

    Attributes
    ----------
    defaultFieldValues: list[tuple[Field, float]]
    regions : OpenFOAMList
    """
    def __init__(self, region: str=""):
        super().__init__("setFieldsDict", "system", region)

        self.defaultFieldValues: list[tuple[Field, float]] = []
        self.regions: OpenFOAMList = OpenFOAMList(SetFieldRegion, "regions")


    def append(self, item):
        check_type("item", item, SetFieldRegion)
        self.regions.append(item)

    @property
    def is_empty(self):
        return(self.regions.is_empty)


    def add_default_values(self, field: Field, value: float | Vector):
        check_type("field", field, Field)
        check_type("value", value, (int, float, Vector))
        self.defaultFieldValues.append((field, value))


    def add_cylinder_to_cell(
            self,
            fieldValues: list[tuple[Field, float]],
            radius: float,
            point1: Vector=Vector(0, 0, 0),
            point2: Vector=Vector(0, 0, 1),
        ):
        """
        Set cell values in a cylinder

        Parameters
        ----------
        fieldValues : list[tuple[Field, float]]
            List of field with the corresponding value to assign.
        radius : float
            Radius of the cylinder.
        point1 : Vector
            Center of the first face of the cylinder.
        point2 : Vector
            Center of the second face of the cylinder.
        """
        self.append(CylinderToCell(fieldValues, radius, point1, point2))


    def add_cylinder_annulus_to_cell(
            self,
            fieldValues: list[tuple[Field, float]],
            innerRadius: float,
            outerRadius: float,
            point1: Vector=Vector(0, 0, 0),
            point2: Vector=Vector(0, 0, 1),
        ):
        """
        Set cell values in a hollow cylinder

        Parameters
        ----------
        fieldValues : list[tuple[Field, float]]
            List of field with the corresponding value to assign.
        innerRadius : float
            Inner radius of the cylinder.
        outerRadius : float
            Outer radius of the cylinder.
        point1 : Vector
            Center of the first face of the cylinder.
        point2 : Vector
            Center of the second face of the cylinder.
        """
        self.append(CylinderAnnulusToCell(
            fieldValues, innerRadius, outerRadius, point1, point2
        ))


    def add_sphere_to_cell(
            self,
            fieldValues: list[tuple[Field, float]],
            outerRadius: float,
            innerRadius: float=0,
            origin: Vector=Vector(0, 0, 0),
        ):
        """
        Set cell values in a hollow sphere

        Parameters
        ----------
        fieldValues : list[tuple[Field, float]]
            List of field with the corresponding value to assign.
        outerRadius : float
            Outer radius of the sphere.
        innerRadius : float
            Inner radius of the sphere (default `0`).
        origin : Vector
            Center of the sphere (default center on the
            origin `Vector(0, 0, 0)`).
        """
        self.append(SphereToCell(
            fieldValues, innerRadius, outerRadius, origin
        ))


    def add_box_to_cell(
            self,
            fieldValues: list[tuple[Field, float]],
            lowCorner: Vector,
            highCorner: Vector,
        ):
        """
        Set cell values in a cylinder

        Parameters
        ----------
        fieldValues : list[tuple[Field, float]]
            List of field with the corresponding value to assign.
        lowCorner : Vector
            Coordinates of the lower corner of the parallelepiped.
        highCorner : Vector
            Coordinates of the upper corner of the parallelepiped.
        """
        self.append(BoxToCell(fieldValues, lowCorner, highCorner))


    def add_zone_to_cell(
            self,
            fieldValues: list[tuple[Field, float]],
            zoneName: str
        ):
        """
        Set cell values in a zone.

        Parameters
        ----------
        fieldValues : list[tuple[Field, float]]
            List of field with the corresponding value to assign.
        zoneName : str
            Name of the zone.
        """
        self.append(ZoneToCell(fieldValues, zoneName))


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter(
            "defaultFieldValues", getFieldValuesAsText(self.defaultFieldValues),
            isAddExtraLine=True
        )

        text += f"{self.regions}"

        return(text)
