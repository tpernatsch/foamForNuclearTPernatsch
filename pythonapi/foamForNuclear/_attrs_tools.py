from __future__ import annotations

from typing import get_origin, get_args, Union, Dict
from types import UnionType
import attrs as attr
from attrs import validators as v
from functools import partial
from foamForNuclear.common import List, Vector, _resolve_type, CheckedOpenFOAMDict, CheckedOpenFOAMList
from collections.abc import Iterable


# CheckedList and CheckedDict assumed imported / defined above

try:
    from typing import dataclass_transform  # Python 3.11+
except ImportError:
    from typing_extensions import dataclass_transform  # pip install typing_extensions


def _make_list_converter(elem_type, container_type):
    """
    Generic converter that:
      - accepts None, container_type, list, tuple, etc.
      - always returns container_type([...]) internally
      - optionally enforces elem_type on entries
    """
    def _conv(value):
        if value is None:
            return container_type()

        # Already the right container
        if isinstance(value, container_type):
            seq = list(value)
        # Plain list/tuple/iterable
        elif isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
            seq = list(value)
        # Single scalar -> treat as [scalar]
        else:
            seq = [value]

        # Type check elements if elem_type is concrete
        if isinstance(elem_type, type):
            for x in seq:
                if not isinstance(x, elem_type):
                    raise TypeError(
                        f"Expected {elem_type.__name__} in list, got {type(x).__name__}"
                    )

        return container_type(seq)

    return _conv


def auto_type_validator(_cls, fields):
    """
    Field transformer that auto-attaches a type validator to every field:
      - str, int, float, bool                -> instance_of(...)
      - Optional[T] / T | None (PEP 604)     -> optional(instance_of(T))
      - List[T] / list[T]                    -> CheckedOpenFOAMList[T]
      - Dict[K, V] / dict[K, V]              -> CheckedOpenFOAMDict[K, V]
      - Any other concrete class             -> instance_of(class)
    Skips fields that already have a validator.
    """
    out = []
    for f in fields:
        # Respect explicit validators
        if f.validator is not None:
            out.append(f)
            continue

        t = _resolve_type(f.type, _cls)
        if t is None:
            out.append(f)
            continue

        origin = get_origin(t)
        args   = get_args(t)

        # --- SPECIAL CASE: Dict[K, V] -> CheckedOpenFOAMDict ----------------
        if origin in (dict, Dict):
            key_type   = args[0] if len(args) >= 1 else object
            value_type = args[1] if len(args) >= 2 else object

            key_validator = (
                v.instance_of(key_type) if isinstance(key_type, type)
                else (lambda *_: None)
            )
            value_validator = (
                v.instance_of(value_type) if isinstance(value_type, type)
                else (lambda *_: None)
            )

            mapping_validator = v.instance_of(CheckedOpenFOAMDict)

            validator = v.deep_mapping(
                key_validator=key_validator,
                value_validator=value_validator,
                mapping_validator=mapping_validator,
            )

            def _make_checkeddict_converter(k_t, v_t, field_name: str):
                def conv(value):
                    if isinstance(value, CheckedOpenFOAMDict):
                        return value
                    if value is None:
                        return CheckedOpenFOAMDict(k_t, v_t, field_name)
                    if isinstance(value, dict):
                        return CheckedOpenFOAMDict(k_t, v_t, field_name, value)
                    try:
                        return CheckedOpenFOAMDict(k_t, v_t, field_name, dict(value))
                    except TypeError:
                        raise TypeError(
                            f"Cannot convert value {value!r} to "
                            f"CheckedOpenFOAMDict[{k_t}, {v_t}] for field '{field_name}'"
                        )
                return conv

            converter = _make_checkeddict_converter(key_type, value_type, f.name)

            f = f.evolve(
                validator=validator,
                converter=converter if f.converter is None else f.converter,
            )
            out.append(f)
            continue
        # -----------------------------------------------------------------

        # --- SPECIAL CASE: list[T] | List[T]  -> CheckedOpenFOAMList --------
        if origin in (list, List):
            elem_type = args[0] if args else object

            member_validator = (
                v.instance_of(elem_type) if isinstance(elem_type, type)
                else (lambda *_: None)
            )

            iterable_validator = v.instance_of(CheckedOpenFOAMList)

            validator = v.deep_iterable(
                member_validator=member_validator,
                iterable_validator=iterable_validator,
            )

            def _make_checkedlist_converter(elem_t, field_name: str):
                def conv(value):
                    if isinstance(value, CheckedOpenFOAMList):
                        return value
                    if value is None:
                        return CheckedOpenFOAMList(elem_t, field_name)
                    if isinstance(value, (list, tuple, set)):
                        return CheckedOpenFOAMList(elem_t, field_name, value)
                    if isinstance(value, elem_t):
                        return CheckedOpenFOAMList(elem_t, field_name, [value])
                    try:
                        return CheckedOpenFOAMList(elem_t, field_name, list(value))
                    except TypeError:
                        raise TypeError(
                            f"Cannot convert value {value!r} to "
                            f"CheckedOpenFOAMList[{elem_t}] for field '{field_name}'"
                        )
                return conv

            converter = _make_checkedlist_converter(elem_type, f.name)

            f = f.evolve(
                validator=validator,
                converter=converter if f.converter is None else f.converter,
            )
            out.append(f)
            continue
        # -----------------------------------------------------------------

        # Basic scalar types
        if t in (str, int, float, bool):
            f = f.evolve(validator=v.instance_of(t))

        # Optional[...] / Union[..., None]
        elif origin in (Union, UnionType):
            has_None = False
            non_none_args = []
            for a in args:
                if a is type(None):
                    has_None = True
                else:
                    non_none_args.append(a)

            # ---- Special case: Union includes dict[K,V] (e.g. dict[str,int] | int, or dict[...] | None) ----
            dict_arg = None
            other_args = []
            for a in non_none_args:
                if get_origin(a) in (dict, Dict):
                    dict_arg = a
                else:
                    other_args.append(a)

            if dict_arg is not None:
                d_args = get_args(dict_arg)
                key_type   = d_args[0] if len(d_args) >= 1 else object
                value_type = d_args[1] if len(d_args) >= 2 else object

                key_validator = (
                    v.instance_of(key_type) if isinstance(key_type, type)
                    else (lambda *_: None)
                )
                value_validator = (
                    v.instance_of(value_type) if isinstance(value_type, type)
                    else (lambda *_: None)
                )
                mapping_validator = v.instance_of(CheckedOpenFOAMDict)

                dict_validator = v.deep_mapping(
                    key_validator=key_validator,
                    value_validator=value_validator,
                    mapping_validator=mapping_validator,
                )

                other_runtime_types = [a for a in other_args if isinstance(a, type)]
                other_validator = None
                if other_runtime_types:
                    other_validator = v.instance_of(
                        tuple(other_runtime_types) if len(other_runtime_types) > 1 else other_runtime_types[0]
                    )

                if other_validator is not None:
                    base_validator = v.or_(dict_validator, other_validator)
                else:
                    base_validator = dict_validator

                validator = v.optional(base_validator) if has_None else base_validator

                field_name = f.name  # bind now (avoid late-binding closure bug)

                def conv(value):
                    if value is None:
                        return None if has_None else value
                    if isinstance(value, CheckedOpenFOAMDict):
                        return value
                    if isinstance(value, dict):
                        return CheckedOpenFOAMDict(key_type, value_type, field_name, value)
                    try:
                        as_dict = dict(value)
                    except TypeError:
                        return value  # e.g. int -> keep as-is
                    else:
                        return CheckedOpenFOAMDict(key_type, value_type, field_name, as_dict)

                f = f.evolve(
                    validator=validator,
                    converter=conv if f.converter is None else f.converter,
                )
                out.append(f)
                continue
            # --------------------------------------------------------------------------------------------

            # ---- Special case: Union includes list[T] (e.g. list[str] | int, or list[...] | None) ----
            list_arg = None
            other_args = []
            for a in non_none_args:
                if get_origin(a) in (list, List):
                    list_arg = a
                else:
                    other_args.append(a)

            if list_arg is not None:
                l_args = get_args(list_arg)
                elem_type = l_args[0] if l_args else object

                member_validator = (
                    v.instance_of(elem_type) if isinstance(elem_type, type)
                    else (lambda *_: None)
                )
                iterable_validator = v.instance_of(CheckedOpenFOAMList)

                list_validator = v.deep_iterable(
                    member_validator=member_validator,
                    iterable_validator=iterable_validator,
                )

                other_runtime_types = [a for a in other_args if isinstance(a, type)]
                other_validator = None
                if other_runtime_types:
                    other_validator = v.instance_of(
                        tuple(other_runtime_types) if len(other_runtime_types) > 1 else other_runtime_types[0]
                    )

                if other_validator is not None:
                    base_validator = v.or_(list_validator, other_validator)
                else:
                    base_validator = list_validator

                validator = v.optional(base_validator) if has_None else base_validator

                field_name = f.name  # bind now (avoid late-binding closure bug)

                def conv(value):
                    if value is None:
                        return None if has_None else value
                    if isinstance(value, CheckedOpenFOAMList):
                        return value
                    if isinstance(value, (list, tuple, set)):
                        return CheckedOpenFOAMList(elem_type, field_name, value)
                    if isinstance(elem_type, type) and isinstance(value, elem_type):
                        return CheckedOpenFOAMList(elem_type, field_name, [value])
                    try:
                        return CheckedOpenFOAMList(elem_type, field_name, list(value))
                    except TypeError:
                        return value  # e.g. int -> keep as-is

                f = f.evolve(
                    validator=validator,
                    converter=conv if f.converter is None else f.converter,
                )
                out.append(f)
                continue
            # --------------------------------------------------------------------------------------------

            # ---- Special case: Optional[list[T]] / list[T] | None ----
            if has_None and len(non_none_args) == 1:
                a0 = non_none_args[0]
                a0_origin = get_origin(a0)
                a0_args   = get_args(a0)

                if a0_origin in (list, List):
                    elem_type = a0_args[0] if a0_args else object

                    member_validator = (
                        v.instance_of(elem_type) if isinstance(elem_type, type)
                        else (lambda *_: None)
                    )
                    iterable_validator = v.instance_of(CheckedOpenFOAMList)

                    base_validator = v.deep_iterable(
                        member_validator=member_validator,
                        iterable_validator=iterable_validator,
                    )

                    validator = v.optional(base_validator)

                    def conv(value):
                        if value is None:
                            return None
                        if isinstance(value, CheckedOpenFOAMList):
                            return value
                        if isinstance(value, (list, tuple, set)):
                            return CheckedOpenFOAMList(elem_type, f.name, value)
                        if isinstance(elem_type, type) and isinstance(value, elem_type):
                            return CheckedOpenFOAMList(elem_type, f.name, [value])
                        try:
                            return CheckedOpenFOAMList(elem_type, f.name, list(value))
                        except TypeError:
                            raise TypeError(
                                f"Cannot convert value {value!r} to "
                                f"CheckedOpenFOAMList[{elem_type}] for field '{f.name}'"
                            )

                    f = f.evolve(
                        validator=validator,
                        converter=conv if f.converter is None else f.converter,
                    )
                    out.append(f)
                    continue

                # ---- Special case: Optional[dict[K,V]] / dict[K,V] | None ----
                if a0_origin in (dict, Dict):
                    key_type   = a0_args[0] if len(a0_args) >= 1 else object
                    value_type = a0_args[1] if len(a0_args) >= 2 else object

                    key_validator = (
                        v.instance_of(key_type) if isinstance(key_type, type)
                        else (lambda *_: None)
                    )
                    value_validator = (
                        v.instance_of(value_type) if isinstance(value_type, type)
                        else (lambda *_: None)
                    )
                    mapping_validator = v.instance_of(CheckedOpenFOAMDict)

                    base_validator = v.deep_mapping(
                        key_validator=key_validator,
                        value_validator=value_validator,
                        mapping_validator=mapping_validator,
                    )

                    validator = v.optional(base_validator)

                    def conv(value):
                        if value is None:
                            return None
                        if isinstance(value, CheckedOpenFOAMDict):
                            return value
                        if isinstance(value, dict):
                            return CheckedOpenFOAMDict(key_type, value_type, f.name, value)
                        try:
                            return CheckedOpenFOAMDict(key_type, value_type, f.name, dict(value))
                        except TypeError:
                            raise TypeError(
                                f"Cannot convert value {value!r} to "
                                f"CheckedOpenFOAMDict[{key_type}, {value_type}] for field '{f.name}'"
                            )

                    f = f.evolve(
                        validator=validator,
                        converter=conv if f.converter is None else f.converter,
                    )
                    out.append(f)
                    continue

            # ---- Fallback: Union of plain runtime types (int|float|None etc.) ----
            members = []
            for a in non_none_args:
                if isinstance(a, type):
                    members.append(a)

            if members:
                base = v.instance_of(tuple(members) if len(members) > 1 else members[0])
                f = f.evolve(validator=v.optional(base) if has_None else base)
            else:
                print(f"WARNING: No validation on field {f.name}")

        # Any other simple concrete type
        elif isinstance(t, type):
            f = f.evolve(validator=v.instance_of(t))

        else:
            print(f"WARNING: No validation on field {f.name}")

        out.append(f)

    return out


def _validate_vec3(instance, attr, value):
    if not isinstance(value, (list, List, Vector)):
        raise TypeError(f"{attr.name} must be a list, List or ffn.Vector. Got {type(value).__name__}")
    if len(value) != 3:
        raise ValueError(f"{attr.name} must have length 3, got {len(value)}")
    for i, v in enumerate(value):
        if not isinstance(v, (int, float)):
            raise TypeError(f"{attr.name}[{i}] must be a number, got {type(v).__name__}: {v!r}")

def call_method_on_change(method_name: str, *field_names: str):
    wanted = set(field_names)
    def _hook(self, attribute: attr.Attribute, value):
        if attribute.name in wanted:
            m = getattr(self, method_name, None)
            if m is not None:
                m(value)            # pass NEW value (on_setattr runs before assign)
        return value
    return _hook

def _is_number_when_flag(flag_name: str):
    def _validator(inst, attribute, value):
        if getattr(inst, flag_name):
            if value is None:
                raise ValueError(f"{attribute.name} is required when {flag_name} is True")
            if not isinstance(value, (int, float)):
                raise TypeError(f"{attribute.name} must be int or float when {flag_name} is True")
    return _validator

def _propagate_region_to(inst, value, *child_names: str) -> None:
    """Write `value` to child.region for each child if present."""
    for name in child_names:
        obj = getattr(inst, name, None)
        if obj is not None and hasattr(obj, "region"):
            obj.region = value

from typing_extensions import dataclass_transform  # or from typing in Py3.11+

def mirror_to_dict(inst, attribute, value):    
    # Skip internal attrs fields (like _of itself)
    if attribute.metadata.get("ffn_internal", False):
        return value

    # Skip name (we don't want it in the dict)
    if attribute.name == "name":
        # but keep _of.name synced
        if hasattr(inst, "_of"):
            inst._of.name = value
        return value

    # NOTE: inst._of is an OpenFOAMDict (or wrapper) that supports pop/setitem
    if hasattr(inst, "_of"):
        if value is None:
            inst._of.pop(attribute.name, None)
        else:
            inst._of[attribute.name] = value

    return value

@dataclass_transform(field_specifiers=(attr.field,))
def ffn_define(_cls=None, **kwargs):
    def wrap(cls):
        return attr.define(
            slots=True,
            on_setattr=[attr.setters.convert, attr.setters.validate, mirror_to_dict],
            field_transformer=auto_type_validator, 
            repr=False,
            kw_only=True
        )(cls)
    return wrap if _cls is None else wrap(_cls)

# tiny converters

def _to_Vector(x) -> Vector:
    """Accept Vector, List, list/tuple, or numpy ndarray; return ffn Vector wrapper."""
    if isinstance(x, Vector):
        return x 
    # allow List or sequences
    if isinstance(x, (List, list, tuple)):
        if(len(x)==3.0):
            if all(isinstance(x_i, (float, int)) for x_i in x):
                return Vector(x[0], x[1], x[2])
            else:
                raise TypeError(f"For Vector {x}, expected either Vector or List, list/tuple, or numpy.ndarray of size 3.0")
        else:
            raise TypeError(f"For Vector {x}, expected either Vector or List, list/tuple, or numpy.ndarray of size 3.0")

    # allow plain sequences
    if isinstance(x, (list, tuple)):
        if(len(x)==3.0):
            return Vector(x)
    # # allow numpy arrays without hard-depending on numpy
    # try:
    #     import numpy as np  # local import: optional dependency
    #     if isinstance(x, np.ndarray):
    #         return List(x.tolist())
    # except Exception:
    #     pass
    raise TypeError("Expected either Vector or List, list/tuple, or numpy.ndarray of size 3.0")

def _to_List_any(x) -> List:
    """Accept List, list/tuple, or numpy ndarray; return ffn List wrapper."""
    if isinstance(x, List):
        return x
    # allow plain sequences
    if isinstance(x, (list, tuple)):
        return List(x)
    # allow numpy arrays without hard-depending on numpy
    try:
        import numpy as np  # local import: optional dependency
        if isinstance(x, np.ndarray):
            return List(x.tolist())
    except Exception:
        pass
    raise TypeError("Expected List, list/tuple, or numpy.ndarray")

def _to_List_float(x) -> List:
    L = _to_List_any(x)
    # coerce elements to float (or leave as is if you prefer)
    return List([float(v) for v in L])

def _to_List_str(x) -> List:
    L = _to_List_any(x)
    return List([str(v) for v in L])

def _is_iterable_not_str(x):
    try:
        iter(x)
    except TypeError:
        return False
    return not isinstance(x, (str, bytes))

def _repr_trunc(obj, maxlen=120):
    s = repr(obj)
    return s if len(s) <= maxlen else s[: maxlen - 3] + "..."

def _to_List_tuple_float(x) -> List:
    """
    Convert x to your List wrapper of tuples[float, ...].
    Accepts: your List, list/tuple of iterables, numpy arrays.
    """
    L = _to_List_any(x)  # reuses your existing helper -> returns List or wraps to List
    out = []
    for i, item in enumerate(L):
        if not _is_iterable_not_str(item):
            raise TypeError(f"List \"{_repr_trunc(x)}\" must be a tuple-list! The element [{i}] must be an iterable of numbers (e.g., tuple[float, ...]); "
                f"got \"{item}\" of type {type(item)}!")
        try:
            tup = tuple(float(v) for v in item)
        except (TypeError, ValueError):
            raise TypeError(f"List \"{_repr_trunc(x)}\" must be a tuple-list! The element [{i}] contains non-numeric values: {item!r}")
        out.append(tup)
    return List(out)

def _to_List_tuple2_float(x) -> List:
    """
    Like above but enforces exactly 2 elements per tuple.
    """
    L = _to_List_any(x)
    out = []
    for i, item in enumerate(L):
        if not _is_iterable_not_str(item):
            raise TypeError(f"List \"{_repr_trunc(x)}\" must be a tuple-list! The element [{i}] must be a 2-tuple/list of numbers, got {type(item)!r}")
        try:
            a, b = item  # will raise if len != 2
        except Exception:
            raise TypeError(f"List \"{_repr_trunc(x)}\" must be a tuple-list! The element [{i}] must have length 2, got {item!r}")
        try:
            tup = (float(a), float(b))
        except (TypeError, ValueError):
            raise TypeError(f"List \"{_repr_trunc(x)}\" must be a tuple-list! The element [{i}] contains non-numeric values: {item!r}")
        out.append(tup)
    return List(out)