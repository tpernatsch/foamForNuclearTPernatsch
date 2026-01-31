from copy import copy
import numpy as np
import base64
import errno
import io, requests
import shutil
from PIL import Image as im
import matplotlib.pyplot as plt

from foamForNuclear.checkvalue import CheckedList, check_type


tab = "    "


class Color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


class Keyword:
    ERROR = f"{Color.RED}Error{Color.END}"
    WARNING = f"{Color.YELLOW}Warning{Color.END}"
    DONE = f"{Color.GREEN}Done{Color.END}"


def underline(text: str) -> str:
    return(f"{Color.UNDERLINE}{text}{Color.END}")

def bold(text: str) -> str:
    return(f"{Color.BOLD}{text}{Color.END}")


class OpenFOAMDict(dict):
    """
    Overload of the Python dictionary. Print in OpenFOAM format.
    """

    def __init__(self, items: dict=None, name: str=None):
        if (items is not None):
            for key, item in items.items():
                self[key] = item

        self.name = name


    def __repr__(self, depth: int=0):
        text = ""
        if (self.name is not None):
            text += f"{self.name}"

        text += f"\n{depth*tab}" + "{\n"
        for key, item in self.items():
            key = format_to_openfoam_regex(key)
            # if ("," in key):
            #     key = "\""+key+"\""
            if (type(item) == bool):
                item = "true" if item else "false"

            if (isinstance(item, OpenFOAMDict)):
                text += f"{(depth+1)*tab}{key} {item.__repr__(depth=depth+1)}\n"
            elif (isinstance(item, (OpenFOAMList, OpenFOAMListDict))):
                text += f"{item.__repr__(depth=depth+1)}\n"
            elif (isinstance(item, Table)):
                text += f"{(depth+1)*tab}{key:15} {item.__repr__(depth=depth+1)};\n"
            else:
                if (item is not None or item == ""):
                    text += f"{(depth+1)*tab}{key:15} {item};\n"
                else:
                    text += f"{(depth+1)*tab}{key};\n"

        text += depth*tab + "}\n"
        return(text)

    @property
    def is_empty(self):
        return(len(list(self.keys())) == 0)


class OpenFOAMListDict(CheckedList):
    def __init__(self, expected_type, name, items=None):
        super().__init__(expected_type, name, items)

    def __repr__(self, depth = 0):
        text = (depth*tab) + self.name + "\n" + (depth*tab) + "{\n"
        for zone in self:
            if (isinstance(zone, OpenFOAMListDict)):
                text += f"{depth*tab}{zone.__repr__(depth=depth+1)}\n"
            else:
                text += f"{(depth+1)*tab}{zone.__repr__(depth=depth+1)}\n"
        text += (depth*tab) + "}\n"
        return(text)

    @property
    def is_empty(self):
        return(self.__len__() == 0)


class OpenFOAMList(CheckedList):
    def __init__(self, expected_type, name, items=None):
        super().__init__(expected_type, name, items)

    def __repr__(self, depth = 0):
        text = f"{(depth*tab)}{self.name}\n{depth*tab}(\n"
        for zone in self:
            if (isinstance(zone, str)):
                text += f"{(depth+1)*tab}{zone}\n"
            else:
                text += f"{(depth+1)*tab}{zone.__repr__(depth=depth+1)}\n"
        text += (depth*tab) + ");\n"
        return(text)

    @property
    def is_empty(self):
        return(self.__len__() == 0)


class Vector:
    """
    Vector of dimension 3

    Parameters
    ----------
    x : float
    y : float
    z : float
    """

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return(f"({self.x:g} {self.y:g} {self.z:g})")

    def __eq__(self, value):
        eps = 1e-6
        return(
            abs(self.x - value.x) <= eps and
            abs(self.y - value.y) <= eps and
            abs(self.z - value.z) <= eps
        )

    def __hash__(self):
        """
        To be used with care
        """
        eps = 1e-8
        return hash((
            round(self.x / eps),
            round(self.y / eps),
            round(self.z / eps),
        ))

    def __add__(self, rhs):
        return(Vector(
            x=self.x + rhs.x,
            y=self.y + rhs.y,
            z=self.z + rhs.z
        ))

    def __sub__(self, rhs):
        return(Vector(
            x=self.x - rhs.x,
            y=self.y - rhs.y,
            z=self.z - rhs.z
        ))

    def __neg__(self):
        return(Vector(
            x=-self.x,
            y=-self.y,
            z=-self.z
        ))

    def dot(self, rhs) -> float:
        """
        Dot product
        """
        return(self.x * rhs.x + self.y * rhs.y + self.z * rhs.z)

    def cross(self, rhs):
        """
        Cross product
        """
        return(Vector(
            x=self.y*rhs.z - self.z*rhs.y,
            y=self.z*rhs.x - self.x*rhs.z,
            z=self.x*rhs.y - self.y*rhs.x
        ))

    def norm(self) -> float:
        return(np.sqrt(self.x**2 + self.y**2 + self.z**2))

    def normalize(self, norm: float=1) -> None:
        if (norm != 0):
            norm = self.norm() / norm

            self.x /= norm
            self.y /= norm
            self.z /= norm
        else:
            self.x = 0
            self.y = 0
            self.z = 0

    def flip(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z

    def translate(self, dx: float=0, dy: float=0, dz: float=0) -> None:
        self.x += dx
        self.y += dy
        self.z += dz

    def rotateX(self, theta: float=0) -> None:
        py = self.y
        self.y =  py * np.cos(theta) - self.z * np.sin(theta)
        self.z =  py * np.sin(theta) + self.z * np.cos(theta)

    def rotateY(self, theta: float=0) -> None:
        px = self.x
        self.x =  px * np.cos(theta) - self.z * np.sin(theta)
        self.z =  px * np.sin(theta) + self.z * np.cos(theta)

    def rotateZ(self, theta: float=0) -> None:
        px = self.x
        self.x =  px * np.cos(theta) - self.y * np.sin(theta)
        self.y =  px * np.sin(theta) + self.y * np.cos(theta)


class Tensor:
    def __init__(self):
        pass


class List(list):
    """
    List of elements
    """
    def __init__(self, items):
        super().__init__()

        if items is not None:
            for item in items:
                self.append(item)

    def __repr__(
            self,
            depth: int=0,
            noBreak: bool=False,
            nCols: int=1,
            isAddLength: bool=False
        ):
        textLength = ""
        if (isAddLength):
            textLength = str(len(self))

        if (noBreak or len(self) <= 3):
            if (all([isinstance(item, (float, int)) for item in self])):
                return(f"{textLength} ({' '.join([f'{e:.6g}' for e in self])})")
            else:
                return(f"{textLength} ({' '.join([f'{e}' for e in self])})")


        text = f"\n{depth*tab}{textLength}\n{depth*tab}(\n"
        for i, value in enumerate(self):
            if (i % nCols == 0):
                text += (depth+1)*tab

            if (all([isinstance(item, (float, int)) for item in self])):
                text += f"{value:.6g}"
            else:
                text += f"{value}"

            if (i % nCols == nCols-1):
                text += "\n"
            else:
                text += " "

        text += f"{depth*tab})"
        return(text)

    @property
    def is_empty(self):
        return(self.__len__() == 0)


class NonUniformList(List):
    def __init__(self, items: list[(float | int)]):
        super().__init__(items)

    def __repr__(self, depth = 0, noBreak = False):
        return(f"nonuniform List<scalar> {len(self)} " + super().__repr__(depth=depth, noBreak=noBreak))


class Matrix(list):
    def __init__(self, items):
        super().__init__()

        if items is not None:
            for item in items:
                self.append(item)

    def __repr__(self, depth = 0):
        nCols = 0
        if (not self.is_empty):
            nCols = len(self[0])
        text = f"{len(self)} {nCols} (\n"
        for line in self:
            text += f"{(depth+1)*tab}({' '.join([str(e) for e in line])})\n"
        text += f"{depth*tab})"
        return(text)

    @property
    def is_empty(self):
        return(self.__len__() == 0)



class Polynome(List):
    """
    p(x) = sum(coeffs[i]*x^i)

    Example of use:

        As a list: Polynome([1, 2]) = 1 + 2*x

        As arguments: Polynome(1, 2, -3) = 1 + 2*x - 3*x^2

        As coefficients: Polynome(x0=1, x2=2, x5=5) = 1 + 2*x^2 + 5*x^5

    Parameters
    ----------
    x0 : float
        Coefficient of x^0
    x1 : float
        Coefficient of x^1
    x2 : float
        Coefficient of x^2
    x3 : float
        Coefficient of x^3
    x4 : float
        Coefficient of x^4
    x5 : float
        Coefficient of x^5
    x6 : float
        Coefficient of x^6
    x7 : float
        Coefficient of x^7
    """

    def __init__(self, *args, **kwargs):

        if (len(args) > 0 and isinstance(args[0], (List, list))):
            n = len(args[0])
            items = list(copy(args[0]))
            for _ in range(8-n):
                items.append(0)

        elif (len(args) > 0 and isinstance(args[0], (float, int))):
            n = len(args)
            items = list(copy(args))
            for _ in range(8-n):
                items.append(0)

        elif (len(kwargs) > 0):
            items = [0, 0, 0, 0, 0, 0, 0, 0]
            for key, item in kwargs.items():
                idx = int(key[1])
                items[idx] = item

        super().__init__(items)


    def __repr__(self, depth=0):
        return super().__repr__(depth, noBreak=True)


    def value(self, x):
        return(sum([coeff * pow(x, i) for i, coeff in enumerate(self)]))


class Table(list):
    def __init__(self, items, type: str='table'):
        super().__init__()

        if items is not None:
            for item in items:
                self.append(item)

        self.type = type

    def __repr__(self, depth = 0):
        text = f"{self.type}\n{depth*tab}(\n"
        for item in self:
            if (isinstance(item, (tuple, list))):
                text += f"{(depth+1)*tab}({' '.join([f'{e!r:15}' if isinstance(e, Vector) else f'{e:15}' for e in item])})\n"

        text += f"{depth*tab})"
        return text

    def append(self, item: tuple | list):
        check_type("item", item, (tuple, list))
        check_type("item[0]", item[0], (float, int))
        check_type("item[1]", item[1], (float, int, Vector))
        super().append(item)

    def value(self, t: float):
        # First point in the table
        ti, valuei = self[0]
        if (t <= ti):
            return(valuei)

        # Last point in the table
        tf, valuef = self[-1]
        if (tf <= t):
            return(valuef)

        # Interpolate in the table
        for tuplef in self[1:]:
            tf, valuef = tuplef
            if (ti <= t and t <= tf):
                m = (valuef - valuei) / (tf - ti)
                return(m * (t - ti) + valuei)

            ti, valuef = tuplef

        return(None)

    @property
    def is_empty(self):
        return(self.__len__() == 0)


class Box:
    def __init__(
            self,
            lowCorner: Vector,
            highCorner: Vector
        ):
        self.lowCorner = lowCorner
        self.highCorner = highCorner

    def __repr__(self):
        return(f"{self.lowCorner} {self.highCorner}")


def convertTimeUnit(time: float, isShortUnit: bool=True) -> str:
    """
    Convert time in seconds to a more readable unit
    """
    if (time / (365*24*3600) >= 1):
        return(f"{time / (365*24*3600):.3g} {'yr' if isShortUnit else 'years'}")
    elif (time / (30*24*3600) >= 1):
        return(f"{time / (30*24*3600):.3g} {'mo' if isShortUnit else 'months'}")
    elif (time / (7*24*3600) >= 1):
        return(f"{time / (7*24*3600):.3g} {'wk' if isShortUnit else 'weeks'}")
    elif (time / (24*3600) >= 1):
        return(f"{time / (24*3600):.3g} {'d' if isShortUnit else 'days'}")
    elif (time / (3600) >= 1):
        return(f"{time / (3600):.3g} {'h' if isShortUnit else 'hours'}")
    elif (time / (60) >= 1):
        return(f"{time / (60):.3g} {'min' if isShortUnit else 'minutes'}")
    elif (time >= 1):
        return(f"{time:.3g} {'s' if isShortUnit else 'seconds'}")
    elif (time*1e3 >= 1):
        return(f"{time*1e3:.3g} {'ms' if isShortUnit else 'milliseconds'}")
    elif (time*1e6 >= 1):
        return(f"{time*1e6:.3g} {'us' if isShortUnit else 'microseconds'}")

    return(f"{time*1e9:.3g} {'ns' if isShortUnit else 'nanoseconds'}")


def format_to_openfoam_regex(text: str) -> str:
    check_type("text", text, str)

    # Check special OpenFOAM characters
    if (
        "|" in text or
        ".*" in text or
        "," in text
    ):
        # Check quotes are not present
        if (
            "\"" not in text and
            "\'" not in text
        ):
            text = f'"{text}"'

    return(text)


def copyFolder(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno in (errno.ENOTDIR, errno.EINVAL):
            shutil.copy(src, dst)
        else: raise


def generate_mermaid_graph_as_image(graph: str, outputFilename: str) -> bool:
    # Change style
    graph = graph.split("\n")
    graph = graph[:2] + ["config:", "  theme: 'neutral'"] + graph[2:]
    graph = "\n".join(graph)

    # Conversion and Mermaid API request
    graphbytes = graph.encode("utf8")
    base64_bytes = base64.urlsafe_b64encode(graphbytes)
    base64_string = base64_bytes.decode("ascii")
    try:
        img = im.open(io.BytesIO(requests.get('https://mermaid.ink/img/' + base64_string, timeout=10).content))
        plt.imshow(img)
        plt.axis('off') # allow to hide axis
        plt.tight_layout()
        plt.savefig(outputFilename, dpi=600)
        plt.close()
        return(True)
    except requests.exceptions.Timeout:
        print(f"{outputFilename} generation timeout")
        return(False)
    except:
        print("Unexpected error")
        return(False)





def addParameter(paramName, value, indent: int=0, isAddExtraLine: bool=False, none_ok: bool=True):
    if (not none_ok and value is None):
        return("")

    if (type(value) == bool):
        value = 'true' if value else 'false'

    isInclude = paramName == "#include"

    text = f"{indent*tab}{paramName:15} {value}"

    if (not isInclude):
        text += ";"
    text += "\n"

    if (isAddExtraLine):
        text += "\n"
    return(text)


openfoamHeader = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2512                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
"""


def openfoamFileHeader(objectName: str, className: str="dictionary") -> str:
    text  =  "FoamFile\n"
    text +=  "{\n"
    text +=  "    version     2.0;\n"
    text +=  "    format      ascii;\n"
    text += f"    class       {className};\n"
    text += f"    object      {objectName};\n"
    text +=  "}\n"
    text +=  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
    return(text)


openfoamFooterLine = """\n// ************************************************************************* //\n"""