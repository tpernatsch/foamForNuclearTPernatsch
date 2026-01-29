#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import yaml


# ----------------------------
# Small utilities
# ----------------------------

def camel_case(name: str) -> str:
    parts = re.split(r"[^A-Za-z0-9]+", (name or "").strip())
    parts = [p for p in parts if p]
    if not parts:
        return "Unnamed"

    def cap(p: str) -> str:
        if re.fullmatch(r"[A-Z0-9]+", p):
            return p
        if p[:1].isupper():
            return p
        return p[:1].upper() + p[1:]

    return "".join(cap(p) for p in parts)


def strip_model_suffix(name: str) -> str:
    s = (name or "").strip()
    if s in ("gapGasModel", "GapGasModel"):
        return s
    if s.lower().endswith("model") and len(s) > 5 and s.lower() != "gapgasmodel":
        return s[:-5]
    return s


def _dedup_preserve_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def indent_block(text: str, n: int = 4) -> str:
    pad = " " * n
    return "\n".join(pad + line if line.strip() else line for line in text.splitlines())


_NUM_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def py_string_literal(declared_type: str, val: Any) -> str:
    t = declared_type.strip()
    if "|" in t:
        t = t.split("|", 1)[0].strip()

    if val is None:
        return "None"

    if t == "bool":
        if isinstance(val, str):
            if val.lower() in ("true", "yes", "on", "1"):
                return "True"
            if val.lower() in ("false", "no", "off", "0"):
                return "False"
        return "True" if bool(val) else "False"

    if t == "int":
        if isinstance(val, str):
            s = val.replace("_", "")
            if not _NUM_RE.match(s) or "." in s or "e" in s.lower():
                raise ValueError(f"Invalid int literal: {val!r}")
            return s
        return str(int(val))

    if t == "float":
        if isinstance(val, str):
            s = val.replace("_", "")
            if not _NUM_RE.match(s):
                raise ValueError(f"Invalid float literal: {val!r}")
            return s
        return repr(float(val))

    if isinstance(val, str):
        s = val.replace("\\", "\\\\").replace("'", "\\'")
        return f"'{s}'"

    return repr(val)


# ----------------------------
# Parsing .H / .C
# ----------------------------

TYPE_NAME_PATTERNS = [
    re.compile(r'TypeName\s*\(\s*"([^"]+)"\s*\)'),
]

ADD_TO_SELECTION_PAT = re.compile(
    r"addToRunTimeSelectionTable\s*\(\s*([A-Za-z_]\w*)\s*,\s*([A-Za-z_]\w*)\s*,"
)

CLASS_NAME_PAT = re.compile(r"\bclass\s+(?P<cls>[A-Za-z_]\w*)\b")

LIST_TYPE_PAT = re.compile(r"^list\s*\[\s*([A-Za-z0-9_]+)\s*\]\s*$")


def _strip_templates(s: str) -> str:
    out: list[str] = []
    depth = 0
    for ch in s:
        if ch == "<":
            depth += 1
            continue
        if ch == ">":
            depth = max(depth - 1, 0)
            continue
        if depth == 0:
            out.append(ch)
    return "".join(out).strip()


def _clean_base_token(tok: str) -> str:
    tok = tok.strip()
    tok = re.sub(r"\b(public|protected|private|virtual)\b", "", tok).strip()
    tok = " ".join(tok.split())
    tok = _strip_templates(tok).strip()
    return tok.strip()


def extract_declared_class_name(h_text: str, stem: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", h_text, flags=re.DOTALL)
    names = [m.group("cls") for m in CLASS_NAME_PAT.finditer(text)]
    if not names:
        return stem

    candidates = [stem, stem[:1].upper() + stem[1:], camel_case(stem)]
    for c in candidates:
        if c in names:
            return c

    return names[0]


def extract_bases_from_h(h_text: str, cpp_class: str) -> list[str]:
    text = re.sub(r"/\*.*?\*/", "", h_text, flags=re.DOTALL)

    for m in CLASS_NAME_PAT.finditer(text):
        if m.group("cls") != cpp_class:
            continue

        tail = text[m.end():]
        end_brace = tail.find("{")
        end_semi = tail.find(";")

        if end_brace == -1 and end_semi == -1:
            head = tail
        else:
            candidates = [i for i in (end_brace, end_semi) if i != -1]
            head = tail[:min(candidates)]

        if ":" not in head:
            return []

        after_colon = head.split(":", 1)[1]
        after_colon = re.sub(r"//.*", "", after_colon)

        parts = [p.strip() for p in after_colon.split(",") if p.strip()]

        bases: list[str] = []
        for p in parts:
            b = _clean_base_token(p)
            if not b:
                continue
            if "::" in b:
                b = b.split("::")[-1]
            bases.append(b)

        return bases

    return []


def extract_type_name_from_h(h_text: str) -> Optional[str]:
    for pat in TYPE_NAME_PATTERNS:
        m = pat.search(h_text)
        if m:
            return m.group(1).strip()
    return None


def extract_selection_base_from_c(c_text: str, cpp_class_name: str) -> Optional[str]:
    for m in ADD_TO_SELECTION_PAT.finditer(c_text):
        base, child = m.group(1), m.group(2)
        if child == cpp_class_name:
            return base
    return None


# ----------------------------
# Model naming rules
# ----------------------------

@dataclass
class CppInfo:
    stem: str
    yaml_path: Path
    header_path: Path
    source_path: Path
    cpp_class_name: str
    type_name: Optional[str]
    bases_from_h: list[str]
    mother_cpp_class: Optional[str]     # mother as STEM (if resolvable)
    selection_base_cpp: Optional[str]   # Base in addToRunTimeSelectionTable(Base, Child,...)


def normalized_type_name(info: CppInfo) -> Optional[str]:
    if not info.type_name:
        return None
    tn = info.type_name.strip()
    if info.selection_base_cpp is None:
        tn = strip_model_suffix(tn)
    return tn


def pick_mother_class(
    bases_from_h: list[str],
    selection_base_cpp: Optional[str],
    cpp_name_to_stem: dict[str, str],
) -> Optional[str]:
    for b in bases_from_h:
        if b in cpp_name_to_stem:
            return cpp_name_to_stem[b]
    if selection_base_cpp and selection_base_cpp in cpp_name_to_stem:
        return cpp_name_to_stem[selection_base_cpp]
    return None


def inheritance_depth(stem: str, mother_map: dict[str, Optional[str]]) -> int:
    d = 0
    seen = {stem}
    cur = stem
    while True:
        m = mother_map.get(cur)
        if not m or m in seen:
            return d
        seen.add(m)
        d += 1
        cur = m


def order_infos_single_inheritance(infos: list[CppInfo], mother_map: dict[str, Optional[str]]) -> list[CppInfo]:
    return sorted(infos, key=lambda i: (inheritance_depth(i.stem, mother_map), i.stem.lower()))


def add_mother_closure(
    infos: list[CppInfo],
    cpp_infos: dict[str, CppInfo],
    mother_map: dict[str, Optional[str]],
) -> list[CppInfo]:
    by_stem: dict[str, CppInfo] = {i.stem: i for i in infos}
    queue = [i.stem for i in infos]
    seen = set(queue)

    while queue:
        s = queue.pop(0)
        m = mother_map.get(s)
        if not m:
            continue
        if m in cpp_infos and m not in seen:
            seen.add(m)
            by_stem[m] = cpp_infos[m]
            queue.append(m)

    original_order = [i.stem for i in infos]
    added = sorted([s for s in by_stem.keys() if s not in set(original_order)], key=str.lower)
    ordered = original_order + added
    return [by_stem[s] for s in ordered]


# ----------------------------
# YAML -> attrs fields
# ----------------------------

PRIMITIVE_TYPE_MAP = {
    "str": "str",
    "string": "str",
    "word": "str",
    "float": "float",
    "double": "float | int",
    "scalar": "float | int",
    "int": "int",
    "label": "int",
    "bool": "bool",
    "boolean": "bool",
    "dictionary": "dict[str, Any]",
    "dict": "dict[str, Any]",
    "table": "Any",
    "vector": "list[int | float]",
    "list": "list[Any]",
}


def yaml_type_to_py(ytype: Optional[str]) -> str:
    if not ytype:
        return "Any"
    t = ytype.strip()
    m = LIST_TYPE_PAT.match(t)
    if m:
        inner = m.group(1)
        inner_py = PRIMITIVE_TYPE_MAP.get(inner, camel_case(inner))
        return f"list[{inner_py}]"
    return PRIMITIVE_TYPE_MAP.get(t, camel_case(t))


def option_to_field_line(opt: dict) -> str:
    key = opt["key"]

    py = (opt.get("py") or {}) if isinstance(opt.get("py"), dict) else {}
    py_type_override = py.get("type")
    py_factory = py.get("factory")

    py_type = str(py_type_override).strip() if py_type_override else yaml_type_to_py(opt.get("type"))

    required = bool(opt.get("required", False))
    has_default = "default" in opt
    default = opt.get("default", None)

    if py_factory:
        return f"{key}: {py_type} = field(factory={py_factory})"

    if required and not has_default:
        return f"{key}: {py_type}"

    if not required and not has_default:
        return f"{key}: {py_type} | None = None"

    if default is None:
        return f"{key}: {py_type} | None = None"

    if isinstance(default, list):
        return (
            f"{key}: {py_type} = field(factory=lambda: {repr(default)})"
            if default
            else f"{key}: {py_type} = field(factory=list)"
        )
    if isinstance(default, dict):
        return (
            f"{key}: {py_type} = field(factory=lambda: {repr(default)})"
            if default
            else f"{key}: {py_type} = field(factory=dict)"
        )
    if isinstance(default, set):
        return (
            f"{key}: {py_type} = field(factory=lambda: {repr(default)})"
            if default
            else f"{key}: {py_type} = field(factory=set)"
        )

    if isinstance(default, (str, int, float, bool)):
        return f"{key}: {py_type} = {py_string_literal(py_type, default)}"

    if py_type.startswith("list[") or py_type == "list[Any]":
        return f"{key}: {py_type} = field(factory=list)"
    if py_type.startswith("dict[") or py_type == "dict[str, Any]":
        return f"{key}: {py_type} = field(factory=dict)"

    return f"{key}: {py_type} | None = None"


def collect_py_imports_from_yaml(y_doc: dict) -> list[str]:
    imps: list[str] = []

    top = y_doc.get("py_imports")
    if isinstance(top, list):
        for s in top:
            if isinstance(s, str) and s.strip():
                imps.append(s.strip())

    for opt in (y_doc.get("options") or []):
        if not isinstance(opt, dict):
            continue
        py = opt.get("py")
        if not isinstance(py, dict):
            continue
        imp_list = py.get("imports")
        if isinstance(imp_list, list):
            for s in imp_list:
                if isinstance(s, str) and s.strip():
                    imps.append(s.strip())

    return _dedup_preserve_order(imps)


# ----------------------------
# Doc building with inheritance
# ----------------------------

def safe_load_yaml(path: Path) -> dict:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
        text = text.replace("\t", "    ")
        return yaml.safe_load(text) or {}
    except yaml.YAMLError as e:
        print(f"[yaml-error] {path}")
        print(f"  {e}")
        return {}


def options_by_key(opts: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for o in opts:
        k = o.get("key")
        if k:
            out[k] = o
    return out


def linearize_mothers(stem: str, mother_map: dict[str, Optional[str]]) -> list[str]:
    chain: list[str] = []
    seen = {stem}
    cur = mother_map.get(stem)
    while cur:
        if cur in seen:
            break
        chain.append(cur)
        seen.add(cur)
        cur = mother_map.get(cur)
    return chain


def build_options_only_docstring(
    y_local: dict,
    yaml_cache: dict[str, dict],
    mother_chain: list[str],
) -> str:
    # Optional top description (only local, not inherited)
    description = (y_local.get("description") or "").strip()

    local_opts = y_local.get("options") or []
    local_map = options_by_key(local_opts)

    inherited_keys_in_order: list[str] = []
    merged: dict[str, dict] = {}

    for anc in reversed(mother_chain):
        y = yaml_cache.get(anc, {}) or {}
        for o in (y.get("options") or []):
            k = o.get("key")
            if not k:
                continue
            merged[k] = o
            if k not in inherited_keys_in_order:
                inherited_keys_in_order.append(k)

    for o in local_opts:
        k = o.get("key")
        if not k:
            continue
        merged[k] = o

    inherited_final = [k for k in inherited_keys_in_order if k in merged and k not in local_map]
    declared_keys = [o.get("key") for o in local_opts if o.get("key")]

    def fmt_type(o: dict) -> str:
        t = (o.get("type") or "Any").strip()
        return t if t else "Any"

    def fmt_meta(o: dict) -> str:
        req = bool(o.get("required", False))
        parts = [f"required: {req}"]
        if "default" in o:
            parts.insert(0, f"default: {o.get('default')!r}")
        return "; ".join(parts)

    def wrap_lines(s: str, width: int = 88) -> list[str]:
        out: list[str] = []
        for raw in (s or "").splitlines():
            line = raw.strip()
            if not line:
                out.append("")
                continue
            while len(line) > width:
                cut = line.rfind(" ", 0, width)
                if cut <= 0:
                    break
                out.append(line[:cut])
                line = line[cut + 1 :].lstrip()
            out.append(line)
        return out

    def render_section(title: str, keys: list[str]) -> list[str]:
        if not keys:
            return []
        lines: list[str] = [title, "-" * len(title)]
        for k in keys:
            o = merged[k]
            t = fmt_type(o)
            desc = (o.get("description") or "").strip()
            lines.append(f"{k} : {t}")

            desc_lines: list[str] = []
            if desc:
                desc_lines.extend(wrap_lines(desc))
            desc_lines.append(f"({fmt_meta(o)})")

            for dl in desc_lines:
                lines.append(f"    {dl}" if dl.strip() else "")
            lines.append("")
        while lines and lines[-1] == "":
            lines.pop()
        return lines

    doc_lines: list[str] = []
    if description:
        doc_lines.extend(wrap_lines(description))
        doc_lines.append("")
        
    doc_lines.extend(render_section("Options (inherited)", inherited_final))
    if doc_lines:
        doc_lines.append("")
    doc_lines.extend(render_section("Options", declared_keys))

    return "\n".join(doc_lines).rstrip()


# ----------------------------
# Code generation
# ----------------------------

PY_HEADER = """\
# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict
"""


def generate_class_code(py_class: str, base_py: str, type_name: str, y_local: dict, doc: str) -> str:
    opts = y_local.get("options") or []
    field_lines = [option_to_field_line(o) for o in opts]

    type_line = f'    TYPE: ClassVar[str] = {py_string_literal("str", type_name)}'
    body_lines = [type_line] + [f"    {ln}" for ln in field_lines]
    body = "\n".join(body_lines)

    if doc:
        doc_block = f'"""\n{indent_block(doc)}\n    """'
        return f"""
@offbeat_define
class {py_class}({base_py}):
    {doc_block}
{body}
""".lstrip("\n")

    return f"""
@offbeat_define
class {py_class}({base_py}):
{body}
""".lstrip("\n")


def generate_family_base_class_code(py_class: str) -> str:
    return f"""
@offbeat_define
class {py_class}(OffbeatDict):
    TYPE: ClassVar[str] = 'none'
""".lstrip("\n")


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("offbeatlib", type=Path, help="Path to offbeatLib root")
    ap.add_argument("--out", type=Path, default=Path("ffn_api_out"), help="Output package dir (root of destination tree)")
    ap.add_argument("--yaml-glob", default="**/*.yaml", help="Glob for yaml under offbeatLib")
    ap.add_argument(
        "--dst-subdir",
        default="offbeat_lib",
        help="Destination subdir under --out where folder packages are created (e.g. offbeat_lib)",
    )
    args = ap.parse_args()

    root: Path = args.offbeatlib.resolve()
    out: Path = args.out.resolve()
    dst_root: Path = out / args.dst_subdir

    yamls = sorted(root.glob(args.yaml_glob))
    if not yamls:
        raise SystemExit(f"No YAML files found under {root} with glob {args.yaml_glob}")

    tmp_declared: dict[str, str] = {}
    for ypath in yamls:
        stem = ypath.stem
        hpath = ypath.with_suffix(".H")
        cpath = ypath.with_suffix(".C")
        if not hpath.exists() or not cpath.exists():
            continue
        htxt = hpath.read_text(encoding="utf-8", errors="ignore")
        tmp_declared[stem] = extract_declared_class_name(htxt, stem)

    cpp_name_to_stem = {cpp_name: stem for stem, cpp_name in tmp_declared.items()}

    cpp_infos: dict[str, CppInfo] = {}
    for ypath in yamls:
        stem = ypath.stem
        hpath = ypath.with_suffix(".H")
        cpath = ypath.with_suffix(".C")
        if not hpath.exists() or not cpath.exists():
            continue

        htxt = hpath.read_text(encoding="utf-8", errors="ignore")
        ctxt = cpath.read_text(encoding="utf-8", errors="ignore")

        cpp_class_name = tmp_declared.get(stem, stem)
        type_name = extract_type_name_from_h(htxt)
        bases = extract_bases_from_h(htxt, cpp_class_name)

        selection_base_cpp = extract_selection_base_from_c(ctxt, cpp_class_name)
        mother = pick_mother_class(bases, selection_base_cpp, cpp_name_to_stem)

        cpp_infos[stem] = CppInfo(
            stem=stem,
            yaml_path=ypath,
            header_path=hpath,
            source_path=cpath,
            cpp_class_name=cpp_class_name,
            type_name=type_name,
            bases_from_h=bases,
            mother_cpp_class=mother,
            selection_base_cpp=selection_base_cpp,
        )

    if not cpp_infos:
        raise SystemExit("No YAML had matching .H and .C files next to it.")

    yaml_cache: dict[str, dict] = {}
    for stem, info in cpp_infos.items():
        yaml_cache[stem] = safe_load_yaml(info.yaml_path)

    mother_map: dict[str, Optional[str]] = {s: info.mother_cpp_class for s, info in cpp_infos.items()}

    families: dict[str, list[CppInfo]] = {}
    for stem, info in cpp_infos.items():
        if not info.selection_base_cpp:
            continue

        fam_key = info.selection_base_cpp
        families.setdefault(fam_key, []).append(info)

        cur = stem
        seen = {cur}
        while True:
            m = mother_map.get(cur)
            if not m or m in seen:
                break
            if m in cpp_infos:
                families[fam_key].append(cpp_infos[m])
            seen.add(m)
            cur = m

    if not families:
        raise SystemExit("No runtime-selectable (addToRunTimeSelectionTable) classes found among YAML docs.")

    for k, infos in list(families.items()):
        uniq: dict[str, CppInfo] = {}
        for i in infos:
            if i.stem not in uniq:
                uniq[i.stem] = i
        families[k] = add_mother_closure(list(uniq.values()), cpp_infos, mother_map)

    cpp_to_py: dict[str, str] = {}
    for stem, info in cpp_infos.items():
        tn = normalized_type_name(info)
        if tn and tn.lower() != "none":
            cpp_to_py[stem] = camel_case(tn)
        else:
            cpp_to_py[stem] = camel_case(strip_model_suffix(stem))

    folder_to_families: dict[Path, list[tuple[str, list[CppInfo]]]] = {}
    for fam_key, infos in families.items():
        base_stem: Optional[str] = cpp_name_to_stem.get(fam_key)
        base_info: Optional[CppInfo] = cpp_infos.get(base_stem) if base_stem else None
        rep = base_info if base_info is not None else infos[0]
        folder_rel = rep.yaml_path.relative_to(root).parent
        folder_to_families.setdefault(folder_rel, []).append((fam_key, infos))

    for folder_rel, fam_list in sorted(folder_to_families.items(), key=lambda kv: str(kv[0]).lower()):
        folder_dir = dst_root / folder_rel
        folder_dir.mkdir(parents=True, exist_ok=True)

        models_file = folder_dir / "models.py"
        init_file = folder_dir / "__init__.py"

        blocks: list[str] = [PY_HEADER.rstrip()]
        defined_names: set[str] = set()
        export_names: list[str] = []

        def emit_alias(alias_name: str, target_name: str) -> None:
            if alias_name == target_name:
                return

            if alias_name not in defined_names:
                blocks.append(f"{alias_name} = {target_name}  # alias")
                defined_names.add(alias_name)

            if alias_name not in export_names:
                export_names.append(alias_name)

            export_names[:] = [n for n in export_names if n != target_name]

        folder_imports: list[str] = []
        for _, infos in fam_list:
            for info in infos:
                y = yaml_cache.get(info.stem, {}) or {}
                folder_imports.extend(collect_py_imports_from_yaml(y))

        uniq_imports = _dedup_preserve_order(folder_imports)
        if uniq_imports:
            blocks.append("\n".join(uniq_imports).rstrip())

        for fam_key, infos in sorted(fam_list, key=lambda x: x[0].lower()):
            fam_base_py = camel_case(fam_key)

            base_stem: Optional[str] = cpp_name_to_stem.get(fam_key)
            base_info: Optional[CppInfo] = cpp_infos.get(base_stem) if base_stem else None

            public_base = camel_case(strip_model_suffix(fam_base_py))

            if base_info is not None:
                y_local = yaml_cache.get(base_info.stem, {}) or {}
                chain = linearize_mothers(base_info.stem, mother_map)
                doc = build_options_only_docstring(y_local, yaml_cache, chain)

                tn = normalized_type_name(base_info)
                type_name = (tn or "none").strip()

                if fam_base_py not in defined_names:
                    blocks.append(
                        generate_class_code(
                            py_class=fam_base_py,
                            base_py="OffbeatDict",
                            type_name=type_name,
                            y_local=y_local,
                            doc=doc,
                        ).rstrip()
                    )
                    defined_names.add(fam_base_py)

                if type_name.lower() != "none":
                    emit_alias(camel_case(type_name), fam_base_py)

            else:
                if fam_base_py not in defined_names:
                    blocks.append(generate_family_base_class_code(fam_base_py).rstrip())
                    defined_names.add(fam_base_py)

            if public_base != fam_base_py:
                emit_alias(public_base, fam_base_py)
            else:
                if fam_base_py not in export_names:
                    export_names.append(fam_base_py)

            for info in order_infos_single_inheritance(infos, mother_map):
                if base_info is not None and info.stem == base_info.stem:
                    continue

                tn = normalized_type_name(info)
                type_name = (tn or "none").strip()

                if type_name.lower() == "none":
                    py_class = camel_case(strip_model_suffix(info.stem))
                else:
                    py_class = camel_case(type_name)

                if py_class in defined_names:
                    continue

                mother_stem = info.mother_cpp_class
                if (not mother_stem) or (mother_stem == info.stem):
                    base_py = fam_base_py
                else:
                    base_py = cpp_to_py.get(mother_stem, fam_base_py)

                chain = linearize_mothers(info.stem, mother_map)
                y_local = yaml_cache.get(info.stem, {}) or {}
                doc = build_options_only_docstring(y_local, yaml_cache, chain)

                blocks.append(generate_class_code(py_class, base_py, type_name, y_local, doc).rstrip())
                defined_names.add(py_class)

                if py_class not in export_names:
                    export_names.append(py_class)

                if info.selection_base_cpp is not None:
                    public_name = camel_case(strip_model_suffix(py_class))
                    if public_name != py_class:
                        emit_alias(public_name, py_class)

        models_file.write_text("\n\n".join(blocks).rstrip() + "\n", encoding="utf-8")

        export_names = _dedup_preserve_order(export_names)

        init_lines: list[str] = [
            "# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.",
            "from __future__ import annotations",
            "",
        ]
        if export_names:
            init_lines.append(f"from .models import {', '.join(export_names)}")
            init_lines.append("")
            init_lines.append(f"__all__ = [{', '.join(repr(n) for n in export_names)}]")
        else:
            init_lines.append("__all__ = []")

        init_file.write_text("\n".join(init_lines).rstrip() + "\n", encoding="utf-8")

        print(f"[gen] {models_file}")
        print(f"[gen] {init_file}  (n={len(export_names)})")

    print(f"\nDone. Generated folder packages in: {dst_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
