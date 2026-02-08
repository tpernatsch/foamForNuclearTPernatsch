from copy import copy
from os import write
import numpy as np
import base64
import errno
import io, requests
import shutil
from PIL import Image as im
import matplotlib.pyplot as plt

import attrs as attr
from attrs import define, field, validators as v
from foamlib import FoamFile
from typing import Mapping, Optional, Union, get_args, get_origin
from types import UnionType

import sys
import os

from foamForNuclear.checkvalue import CheckedDict, CheckedList, check_type


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

    def __init__(self, items: dict=None, name: str | None =None):
        if (items is not None):
            for key, item in items.items():
                self[key] = item

        self.name = name

    def __setitem__(self, key, item):
        if isinstance(item, list) and not isinstance(item, (OpenFOAMList, OpenFOAMListDict, Table)):
            item = List(item)
        super().__setitem__(key, item)

    def export_body_to_foam(self, depth: int = 0) -> str:
        text = ""
        text += f"\n{depth*tab}" + "{\n"
        
        for key, item in self.items():
            key = format_to_openfoam_regex(key)

            if isinstance(item, bool):
                item = "true" if item else "false"

            if isinstance(item, OpenFOAMDict):
                # For the printed name, preference is given to dict.name.
                # If it exists then repr function will write it before the dict body. 
                # If it does not exist, the key (if any) will be used instead
                if item.name is None:
                    text += f"{(depth+1)*tab}{key}{item.__repr__(depth=depth+1)}\n"
                else:
                    text += f"{(depth+1)*tab}{item.__repr__(depth=depth+1)}\n"
            elif isinstance(item, (FoamForNuclearDict, CheckedOpenFOAMDict,)):
                # For the printed name, preference is given to dict.name.
                # If it exists then repr function will write it before the dict body. 
                # If it does not exist, the key (if any) will be used instead
                if item.name is None:
                    text += f"{(depth+1)*tab}{key}{item._as_openfoam_dict().__repr__(depth=depth+1)}\n"
                else:
                    text += f"{(depth+1)*tab}{item._as_openfoam_dict().__repr__(depth=depth+1)}\n"
            elif isinstance(item, (OpenFOAMList, OpenFOAMListDict)):
                text += f"{item.__repr__(depth=depth+1)}\n"
            elif isinstance(item, Table):
                text += f"{(depth+1)*tab}{key:15} {item.__repr__(depth=depth+1)};\n"
            else:
                if item is not None or item == "":
                    text += f"{(depth+1)*tab}{key:15} {item};\n"
                else:
                    text += f"{(depth+1)*tab}{key};\n"

        text += depth*tab + "}\n"
        return text

    def __repr__(self, depth: int = 0) -> str:
        text = ""
        if self.name is not None:
            text += f"{self.name}"
        text += self.export_body_to_foam(depth=depth)
        return text

    @property
    def is_empty(self):
        return(len(list(self.keys())) == 0)
    

class CheckedOpenFOAMDict(CheckedDict):
    """
    CheckedDict that prints in OpenFOAM format using OpenFOAMDict underneath.
    """

    def __init__(self, expected_key_type, expected_value_type,
                 name: str | None = None, items: dict | None = None):
        # reuse CheckedDict logic (types + items + .name)
        super().__init__(expected_key_type, expected_value_type, name, items)

    def export_body_to_foam(self, depth: int = 0) -> str:
        of_dict = OpenFOAMDict(self, name=self.name)
        return of_dict.export_body_to_foam(depth=depth)

    def _as_openfoam_dict(self) -> OpenFOAMDict:
        of_dict = OpenFOAMDict(self, name=self.name)
        return of_dict

    def __repr__(self, depth: int = 0) -> str:
        of_dict = OpenFOAMDict(self, name=self.name)
        return of_dict.__repr__(depth=depth)

    @property
    def is_empty(self):
        return len(self.keys()) == 0
    

# ---------- subclass discovery (no manual registry) ----------
def _all_subclasses(cls):
    for sub in cls.__subclasses__():
        yield sub
        yield from _all_subclasses(sub)

def _find_by_TYPE(base, type_name: Optional[str]):
    want = (type_name or "").strip().lower()
    best = None
    for c in _all_subclasses(base):
        typ = getattr(c, "TYPE", None)
        if isinstance(typ, str) and typ.strip().lower() == want:
            if best is not None:
                raise ValueError(f"Ambiguous type '{type_name}' for {base.__name__}: "
                                 f"{best.__name__} and {c.__name__}")
            best = c
    return best or base

def _resolve_type(t, cls):
    """Resolve string annotations (PEP 563 / future annotations) to real types."""
    if not isinstance(t, str):
        return t
    modns = sys.modules[cls.__module__].__dict__
    try:
        return eval(t, modns, modns)
    except Exception:
        return None
    
#==============================================================================*
# FoamForNuclearDict and other classes

@define
class FoamForNuclearDict:
    """
    Base for attrs-based OFFBEAT/OpenFOAM dictionary sections.

    - Does NOT inherit from dict / OpenFOAMDict (clean completions)
    - Internally builds an OpenFOAMDict for printing/export
    """
    name: str | None = None

    # internal OpenFOAMDict backing store (hidden from user API & attrs-asdict)
    _of: OpenFOAMDict = field(
        factory=OpenFOAMDict,
        repr=False,
        metadata={"ffn_internal": True},
    )

    def __attrs_post_init__(self):
        self._of.name = self.name
        self._update_fields()

    # ----------------------------
    # private helpers
    # ----------------------------
    def _update_fields(self) -> None:
        def _filter(f, v):
            if v is None:
                return False
            if f.metadata.get("ffn_internal", False):
                return False
            if f.name == "name":
                return False
            return True

        # new values from attrs (flat only)
        data = attr.asdict(self, recurse=False, filter=_filter)

        # inject class TYPE as 'type'
        cls_type = getattr(type(self), "TYPE", None)
        if cls_type and "type" not in data:
            data = {"type": cls_type, **data}

        # parent-declared field order (same logic as you had)
        priority: list[str] = []
        for cls in type(self).mro()[1:]:
            anns = getattr(cls, "__annotations__", {})
            for n in anns.keys():
                if n not in priority:
                    priority.append(n)

        items = list(data.items())
        order_in_asdict = {k: i for i, k in enumerate(data.keys())}
        parent_pos = {k: i for i, k in enumerate(priority)}

        def _sort_key(kv):
            k = kv[0]
            if k == "type":
                return (-1, -1)
            if k in parent_pos:
                return (0, parent_pos[k])
            return (1, order_in_asdict[k])

        items.sort(key=_sort_key)

        # rebuild the internal OpenFOAMDict
        self._of.clear()
        self._of.update(items)

    def _as_openfoam_dict(self) -> OpenFOAMDict:
        self._update_fields()
        return self._of

    # dunder is OK (won't clutter completions)
    def __repr__(self, depth: int = 0) -> str:
        return self._as_openfoam_dict().__repr__(depth=depth)

    @property
    def is_empty(self) -> bool:
        # keep your old property name; this is fine as public
        self._update_fields()
        return self._of.is_empty

    # ----------------------------
    # import logic (keep public entry, make core private)
    # ----------------------------
    @classmethod
    def _from_mapping(
        cls,
        section: Mapping[str, object],
        file_path: Optional[str] = None,
        key_path: tuple = (),
    ) -> "FoamForNuclearDict":

        if not bool(section):
            return None

        type_name = section.get("type") if isinstance(section, Mapping) else None
        concrete = _find_by_TYPE(cls, type_name)

        kwargs = {}
        for f in attr.fields(concrete):
            if f.name in ("file_path", "key_path"):
                continue
            if not isinstance(section, Mapping) or f.name not in section:
                continue

            raw = section[f.name]

            t = _resolve_type(f.type, concrete)
            if t is None:
                continue

            origin = get_origin(t)
            args = get_args(t)

            if origin in (Union, UnionType):
                # 1) FoamForNuclearDict takes precedence
                offbeat_cls = next(
                    (a for a in args if isinstance(a, type) and issubclass(a, FoamForNuclearDict)),
                    None,
                )
                if offbeat_cls is not None:
                    kwargs[f.name] = offbeat_cls._from_mapping(
                        raw, file_path=file_path, key_path=key_path + (f.name,)
                    )
                    continue

                # 2) Table second
                if any(a is Table for a in args):
                    if isinstance(raw, list) and raw and isinstance(raw[0], str):
                        kwargs[f.name] = Table(raw, raw[0])
                    else:
                        kwargs[f.name] = Table(raw, "")
                    continue

                # 3) Fallback converters
                for a in args:
                    try:
                        kwargs[f.name] = a(raw)
                        break
                    except Exception:
                        continue
                continue

            if isinstance(t, type) and issubclass(t, FoamForNuclearDict):
                kwargs[f.name] = t._from_mapping(
                    raw, file_path=file_path, key_path=key_path + (f.name,)
                )
                continue

            if t is OpenFOAMDict:
                kwargs[f.name] = OpenFOAMDict(raw)
                continue

            if t is Table:
                kwargs[f.name] = Table(raw)
                continue

            kwargs[f.name] = raw

        return concrete(**kwargs)  # type: ignore[call-arg]

    @classmethod
    def import_from_openfoam(
        cls,
        file_path: str,
        *key_path: str,
    ) -> "FoamForNuclearDict":
        with FoamFile(file_path) as f:
            try:
                section = dict(f[tuple(key_path)])
            except KeyError:
                section = {}

        return cls._from_mapping(section, file_path=file_path, key_path=tuple(key_path))
    

# A helper class to define Sciantix input settings
# Default corresponds to classic LWR settings
@define
class SciantixDict(FoamForNuclearDict):
    """
    iverification: int (0= no verification) 
    igrain_growth: int (1 = ainscough) 
    iinert_gas_behavior: int (1= do IGB) 
    igas_diffusion_coefficient: int (1= Turnbull et al., 1988) 
    iintra_bubble_evolution: int (1=Pizzocri et al., 2018) 
    ibubble_radius: int (1= Olander&Wongy, 2006) 
    iresolution_rate: int (1=Turnbull 1971) 
    itrapping_rate: int (1= Ham 1958) 
    inucleation_rate: int (1= Baker 1971) 
    isolver: int (1= SDA, Pizzocri et al., 2019) 
    iformat_output: int (1 = output.txt, values separated by tabs) 
    igrain_boundary_vacancy_diffusion_coefficient: int (1= Reynolds and Burton, 1979) 
    igrain_boundary_behaviour: int (1= do InterGranularGasBehavior - Pastore et al., 2013; Barani et al., 2017) 
    igrain_boundary_micro_cracking: int (1 = Barani et al., 2017) 
    igrain_recrystallization: int (0 = non active) 
    ifuel_reactor_type: int (0=UO2/PWR) 
    igas_effective_coefficientgas effective: int  (=0) or single atom (=1)
    igas_sweepingadd boundary: int  gas sweeping 
    imicro_cracking_span: int 
    sf_resolution_rate: int 
    sf_trapping_rate: int 
    sf_nucleation_rate: int 
    sf_diffusion_rate: int 
    
    """
    iverification: int = field(default=0, validator=v.instance_of(int))
    igrain_growth: int = field(default=1, validator=v.instance_of(int))
    iinert_gas_behavior: int = field(default=1, validator=v.instance_of(int))
    igas_diffusion_coefficient: int = field(default=1, validator=v.instance_of(int))
    iintra_bubble_evolution: int = field(default=1, validator=v.instance_of(int))
    ibubble_radius: int = field(default=1, validator=v.instance_of(int))
    iresolution_rate: int = field(default=1, validator=v.instance_of(int))
    itrapping_rate: int = field(default=1, validator=v.instance_of(int))
    inucleation_rate:int = field(default=1, validator=v.instance_of(int))
    isolver:int = field(default=1, validator=v.instance_of(int))
    iformat_output:int = field(default=1, validator=v.instance_of(int))
    igrain_boundary_vacancy_diffusion_coefficient:int = field(default=1, validator=v.instance_of(int))
    igrain_boundary_behaviour:int = field(default=1, validator=v.instance_of(int))
    igrain_boundary_micro_cracking:int = field(default=1, validator=v.instance_of(int))
    igrain_recrystallization:int = field(default=0, validator=v.instance_of(int))
    ifuel_reactor_type:int = field(default=0, validator=v.instance_of(int))
    igas_effective_coefficient:int = field(default=1, validator=v.instance_of(int))
    igas_sweeping:int = field(default=1, validator=v.instance_of(int))
    imicro_cracking_span: int = field(default=0, validator=v.instance_of(int))
    sf_resolution_rate: int = field(default=1, validator=v.instance_of(int))
    sf_trapping_rate: int = field(default=1, validator=v.instance_of(int))
    sf_nucleation_rate: int = field(default=1, validator=v.instance_of(int))
    sf_diffusion_rate: int = field(default=1, validator=v.instance_of(int))


class OpenFOAMListDict(CheckedList):
    def __init__(self, expected_type, name, items=None):
        
        super().__init__(expected_type, name, items)

    def __repr__(self, depth = 0):
        text = (depth*tab) + self.name + "\n" + (depth*tab) + "{\n"
        for zone in self:
            if (isinstance(zone, OpenFOAMListDict)):
                text += f"{zone.__repr__(depth=depth+1)}\n"
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


class CheckedOpenFOAMList(CheckedList):
    """
    A CheckedList that prints using the OpenFOAM-style List.__repr__.
    """

    def __init__(self, expected_type, name, items=None):
        # use the same construction / type checking as CheckedList
        super().__init__(expected_type, name, items)

    def __repr__(
        self,
        depth: int = 0,
        noBreak: bool = False,
        nCols: int = 1,
        isAddLength: bool = False,
    ):
        # Reuse the OpenFOAM List formatting on top of our checked content.
        of_list = List(self)
        return of_list.__repr__(
            depth=depth,
            noBreak=noBreak,
            nCols=nCols,
            isAddLength=isAddLength,
        )
    

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
            # Skipt 'table' name if present in items
            if (not(isinstance(items, dict)) and
                    isinstance(items[0],str) and len(items)==2):
                for item in items[1]:
                    self.append(item)
            else:
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

    if (not isInclude and not isinstance(value, dict)):
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

def of_flavour() -> str | None:
    """
    Return "esi" or "foundation" based on WM_PROJECT_VERSION.
    Heuristic:
      - ESI/openfoam.com versions look like: v2312, v2506, ...
      - Foundation versions look like: 10, 11, 12, ...
    """
    ver = os.environ.get("WM_PROJECT_VERSION", "")
    if ver == "":
        return None
    if ver.startswith("v"):
        return "esi"
    return "foundation"