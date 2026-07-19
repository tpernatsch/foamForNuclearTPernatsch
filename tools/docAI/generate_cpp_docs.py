"""
Generate the doc for Sphinx RST files

Now reads per-class YAML docs placed next to the .H/.C files:
  e.g.  myClass.H, myClass.C, myClass.yaml

YAML fields (all optional except type_name is recommended):
  type_name: str
  description: str (can be multiline with | )
  formulation: str (can be multiline with | )
  admonitions:
    - { kind: warning|note|info|tip|..., body: str }
  options:
    - { key: str, type: str, required: bool, default: any, description: str }
  externalOptions:
    - { key: str, type: str, required: bool, default: any, description: str, path?: str }
  usage:
    - { title?: str, comment?: str, snippet: str }

The rest of the script (release notes, index, etc.) is unchanged.
"""

#==============================================================================*
# Imports

import argparse
import os
import re
import sys
from pathlib import Path

# NEW: YAML
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

#==============================================================================*

def write_if_changed(path, text: str) -> None:
    """
    Write text only if it differs from the current file content, so that
    unchanged files keep their mtime and Sphinx skips re-reading them.
    """
    p = Path(path)
    try:
        if p.read_text(encoding="utf-8") == text:
            return
    except OSError:
        pass
    p.write_text(text, encoding="utf-8")

#==============================================================================*

def replaceInlineReference(text):
    """
    Replace "[@REF]" for inline reference
    """
    if ("[@" not in text):
        return(text)

    refValues = []
    for line in text.split("\n"):
        if ("[@" in line):
            lineSplit = line.split("[@")
            for section in lineSplit[1:]:
                refValues.append(section.split("]")[0])

    for refValue in refValues:
        text = text.replace(f"[@{refValue}]", f":ref:`{refValue} <{refValue}>`")

    return(text)


def replaceInlineMath(text):
    """
    Replace "\f$" and "$" for inline math
    """
    # Remove \f$
    mathInlineFlag = True
    while (r"\f$" in text):
        text = re.sub(
            r"\\f\$\s*" if mathInlineFlag else r"\s*\\f\$",
            ":math:`" if mathInlineFlag else "`",
            text,
            1
        )
        mathInlineFlag = not mathInlineFlag

    # Remove $
    refValues = []
    for line in text.split("\n"):
        if ("$" in line):
            lineSplit = line.split("$")
            if (len(lineSplit) % 2 == 0):
                continue
            for section in lineSplit[1::2]:
                refValues.append(section)

    for refValue in refValues:
        # Whitespace inside $...$ would produce ":math:` x `", which is an
        # invalid RST role (rendered as literal text) — strip it.
        if not refValue.strip():
            continue
        text = text.replace(f"${refValue}$", f":math:`{refValue.strip()}`")

    return(text)


def format_tables_for_rst(text: str) -> str:
    tableRow = []
    textSplit = text.split('\n')
    for i, line in enumerate(textSplit):
        if (r"\table" in line):
            nCols = len(textSplit[i+1].split('|'))
            tableRow.append((i+2+len(tableRow), nCols))

    for idx, nCols in tableRow:
        textSplit = textSplit[:idx] + [(nCols-1)*'-|' + '-'] + textSplit[idx:]

    new_text = '\n'.join(textSplit)
    new_text = new_text.replace(r'\table', '')
    new_text = new_text.replace(r'\endtable', '')

    return(new_text)


def format_table_for_rst(text: str) -> str:
    tableRow = []
    textSplit = text.split('\n')
    isRead = False
    new_text = text
    for i, line in enumerate(textSplit):
        if (r"\table" in line):
            nCols = textSplit[i+1].count('|')+1
            textSplit[i] = textSplit[i].replace(r'\table', f"\n.. list-table::\n    :widths:{nCols*' 50'}\n    :header-rows: 1\n", 1)

            isRead = True
            continue

        if (r"\endtable" in line):
            isRead = False

            for idx in tableRow:
                strippedLine = [element.strip() for element in textSplit[idx].strip().split("|")]
                textSplit[idx] = "    * - " + '\n      - '.join(strippedLine)

            textSplit[i] = textSplit[i].replace(r'\endtable', '', 1)
            new_text = '\n'.join(textSplit)

            tableRow = []
            continue

        if (isRead):
            tableRow.append(i)

    return(new_text)


def format_equation_for_rst(text: str) -> str:
    """
    Convert \\f[...\\f] and $$...$$ display math into ``.. math::`` directives.
    Handles equation content on the same line as the delimiters (LLM-generated
    text often writes ``$$ a + b $$`` rather than putting $$ on its own lines).
    """
    def _to_directive(m):
        content = m.group(1).strip()
        if not content:
            return ""
        body = "\n".join(
            f"    {line.strip()}" for line in content.splitlines() if line.strip()
        )
        # Blank lines around the directive keep RST happy mid-paragraph.
        return f"\n\n.. math::\n\n{body}\n"

    text = re.sub(r'\\f\[(.*?)\\f\]', _to_directive, text, flags=re.DOTALL)
    text = re.sub(r'\$\$(.*?)\$\$', _to_directive, text, flags=re.DOTALL)
    return text


_MATH_SEGMENT_RE = re.compile(
    r'(\$\$.*?\$\$|\\f\[.*?\\f\]|\$[^$\n]+\$|:math:`[^`]*`)', re.DOTALL
)


def sanitize_llm_markup(text: str) -> str:
    """
    Deterministic cleanup of common LLM markup mistakes in description and
    formulation fields:
    - \\texttt{...} in prose -> ``...`` (LaTeX commands don't exist in RST);
      inside math -> \\mathtt{...} (MathJax-safe monospace).
    - \\textbf / \\textit in prose -> **...** / *...*
    - missing blank line between a paragraph and a following "- " bullet list
      (without it RST joins the list into the paragraph).
    """
    def _fix(segment: str, in_math: bool) -> str:
        if in_math:
            return segment.replace(r'\texttt', r'\mathtt')
        segment = re.sub(r'\\texttt\{([^{}]*)\}', r'``\1``', segment)
        segment = re.sub(r'\\textbf\{([^{}]*)\}', r'**\1**', segment)
        segment = re.sub(r'\\textit\{([^{}]*)\}', r'*\1*', segment)
        return segment

    parts = _MATH_SEGMENT_RE.split(text)
    text = ''.join(
        _fix(p, _MATH_SEGMENT_RE.fullmatch(p) is not None) for p in parts
    )

    # Drop standalone 'toctreeHere' placeholder lines — leftovers from an old doc
    # workflow that an injector used to replace; they now render as literal text.
    text = '\n'.join(
        ln for ln in text.split('\n')
        if not re.fullmatch(r'\s*toctreehere\s*', ln, re.IGNORECASE)
    )

    lines = text.split('\n')
    fixed: list[str] = []
    for i, ln in enumerate(lines):
        if (ln.lstrip().startswith('- ') and i > 0 and lines[i - 1].strip()
                and not lines[i - 1].lstrip().startswith('- ')):
            fixed.append('')
        fixed.append(ln)
    return '\n'.join(fixed)

# -------------------------------
# OLD .H parsing (now unused)
# (kept for reference, not used)
# -------------------------------

def extract_descriptions_and_usage(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
    except:
        return("", "", "", "")

    description_pattern = re.compile(r"Description\s*\n\s*(.*?)(?=\n\s*(Usage|\\par Options|SourceFiles|Class|\\\*-----|\*/|\\vartable))", re.DOTALL)
    options_pattern = re.compile(r"\\par Options\s*\n\s*(.*?)(?=\n\s*(Usage|SourceFiles|Class|\\\*-----|\*/))", re.DOTALL)
    usage_pattern = re.compile(r"Usage\s*\n\s*(.*?)(?=\n\s*(SourceFiles|Class|\\\*-----|\*/))", re.DOTALL)
    vartable_pattern = re.compile(r"\\vartable.*?\\endvartable", re.DOTALL)

    description = description_pattern.search(content)
    options = options_pattern.search(content)
    usage = usage_pattern.search(content)
    vartable = vartable_pattern.search(content)

    description_text = description.group(1).strip() if description else "No description available."
    options_text = options.group(1).strip() if options else "No options available."
    usage_text = usage.group(1).strip() if usage else "No usage available."
    vartable_text = vartable.group(0).strip() if vartable else "No vartable available."

    return description_text, options_text, usage_text, vartable_text


def is_runTimeSelectable(file_path):
    c_file_path = file_path.replace('.H', '.C')
    if os.path.exists(c_file_path):
        with open(c_file_path, 'r') as file:
            content = file.read()
        return 'addToRunTimeSelectionTable' in content
    return False


def extract_type_name_from_header(header_path: str) -> str | None:
    """
    Extract OpenFOAM runtime TypeName("...") from the header.
    Falls back to defineTypeNameAndDebug(ClassName, ...) if needed.
    """
    try:
        txt = Path(header_path).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None

    m = re.search(r'\bTypeName\s*\(\s*"([^"]+)"\s*\)', txt)
    if m:
        return m.group(1).strip()

    m = re.search(r"\bdefineTypeNameAndDebug\s*\(\s*([A-Za-z_]\w*)\s*,", txt)
    if m:
        return m.group(1).strip()

    return None


# ============================================================
# NEW: YAML-based generation (replaces .H parsing in practice)
# ============================================================

CLASS_NAME_PAT = re.compile(r"\bclass\s+(?P<cls>[A-Za-z_]\w*)\b")


def extract_declared_class_name_from_header(header_path: str) -> str:
    txt = Path(header_path).read_text(encoding="utf-8", errors="ignore")
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.DOTALL)
    m = CLASS_NAME_PAT.search(txt)
    return m.group("cls") if m else Path(header_path).stem


def extract_first_base_from_h(header_path: str, cpp_class: str) -> str | None:
    txt = Path(header_path).read_text(encoding="utf-8", errors="ignore")
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.DOTALL)

    for m in CLASS_NAME_PAT.finditer(txt):
        if m.group("cls") != cpp_class:
            continue

        tail = txt[m.end():]
        end_brace = tail.find("{")
        end_semi = tail.find(";")
        if end_brace == -1 and end_semi == -1:
            head = tail
        else:
            head = tail[:min(i for i in (end_brace, end_semi) if i != -1)]

        if ":" not in head:
            return None

        after = head.split(":", 1)[1]
        after = re.sub(r"//.*", "", after)
        after = re.sub(r"\b(public|protected|private|virtual)\b", "", after)
        after = " ".join(after.split())

        first = after.split(",", 1)[0].strip()
        if "::" in first:
            first = first.split("::")[-1]
        first = re.sub(r"<.*?>", "", first).strip()
        return first or None

    return None


def find_yaml_doc_for_header(header_path: str) -> Path | None:
    hpath = Path(header_path)
    stem = hpath.stem
    candidates = [
        hpath.with_name(f"{stem}.doc.yaml"),
        hpath.with_name(f"{stem}.yaml"),
    ]
    for c in candidates:
        if c.exists() and c.is_file():
            return c
    return None


def _opts_by_key(opts: list[dict]) -> dict[str, dict]:
    out = {}
    for o in opts:
        k = (o.get("key") or "").strip()
        if k:
            out[k] = o
    return out


def build_yaml_cache_and_mothers(ffn_lib_dirs: list[str]):
    yaml_cache = {}
    header_for_class = {}

    for ffn_lib_dir in ffn_lib_dirs:
        src_root = Path(ffn_lib_dir).resolve()

        for h_path in src_root.rglob("*.H"):
            if any(p in {"lnInclude", "Make"} for p in h_path.parts):
                continue

            yaml_path = find_yaml_doc_for_header(str(h_path))
            if yaml_path is None:
                continue

            cpp_class = extract_declared_class_name_from_header(str(h_path))
            yaml_cache[cpp_class] = load_yaml_doc(yaml_path)
            header_for_class[cpp_class] = h_path

    mother_map = {}
    for cls, h_path in header_for_class.items():
        base = extract_first_base_from_h(str(h_path), cls)
        mother_map[cls] = base if (base in yaml_cache) else None

    return yaml_cache, mother_map


def compute_inherited_only_options(
    cls: str,
    yaml_cache: dict[str, dict],
    mother_map: dict[str, str | None],
) -> list[dict]:
    local_opts = yaml_cache.get(cls, {}).get("options") or []
    local_keys = set(_opts_by_key(local_opts).keys())

    merged = {}
    order = []

    cur = mother_map.get(cls)
    seen = {cls}
    chain = []

    while cur and cur not in seen:
        chain.append(cur)
        seen.add(cur)
        cur = mother_map.get(cur)

    for anc in reversed(chain):
        for o in (yaml_cache.get(anc, {}).get("options") or []):
            k = (o.get("key") or "").strip()
            if not k:
                continue
            merged[k] = o
            if k not in order:
                order.append(k)

    inherited_keys = [k for k in order if k not in local_keys]
    return [merged[k] for k in inherited_keys]


def _rst_admonition(kind: str, body: str) -> str:
    """Render a Sphinx admonition."""
    kind = (kind or "note").lower().strip()
    if kind not in {"note", "warning", "tip", "important", "caution", "attention", "hint", "admonition", "danger"}:
        kind = "note"
    lines = [f".. {kind}::", ""]
    for ln in (body or "").rstrip("\n").splitlines():
        lines.append("   " + ln)
    lines.append("")
    return "\n".join(lines)


def _rst_ai_provenance(prov: dict) -> str:
    """
    Render a warning admonition at the end of the page describing AI involvement.
    Returns an empty string if prov is empty or absent.
    """
    if not prov:
        return ""
    model = prov.get("model", "AI")
    date = prov.get("date", "")
    date_str = f" on {date}" if date else ""

    if prov.get("generated"):
        body = (
            f"This page was **fully generated** by ``{model}``{date_str} "
            "and has not been reviewed by a human.\n\n"
            "Equations, options, and descriptions should be verified against the source code."
        )
        return _rst_admonition("warning", body)

    lines = []
    sections = prov.get("modified_sections") or []
    if sections:
        sec_list = ", ".join(f"``{s}``" for s in sections)
        if prov.get("stub_improved"):
            lines.append(f"- **Sections improved from stub:** {sec_list}")
        else:
            lines.append(f"- **Sections modified:** {sec_list}")
    if prov.get("formulation_added"):
        lines.append("- **Formulation section:** added automatically")
    if prov.get("formulation_extended"):
        lines.append("- **Formulation section:** expanded automatically")
    if not lines:
        return ""

    intro = (
        f"Parts of this page were improved from a stub by ``{model}``{date_str}:"
        if prov.get("stub_improved") else
        f"Parts of this page were modified by ``{model}``{date_str}:"
    )
    body = intro + "\n\n" + "\n".join(lines) + "\n\nPlease verify the changes against the source code."
    return _rst_admonition("warning", body)


def _options_need_dict_column(options: list, class_dictionary: str | None) -> bool:
    """True when some option is read from a dictionary file other than the
    class-level default — in that case a per-row Dictionary column is clearer
    than a single 'Set in' note."""
    return any(
        o.get("dictionary") and o.get("dictionary") != class_dictionary
        for o in options
    )


def _rst_options_list_table(options: list, include_path: bool = False,
                            class_dictionary: str | None = None) -> str:
    """Render options (key/[path|location|dictionary]/type/required/default/description) as a list-table."""
    if not options:
        return "No options list available.\n\n"

    # AI-generated yamls may attach a `location` (sub-dictionary) to options;
    # render it as an extra column, reusing the path-column layout.
    include_location = (not include_path) and any(o.get("location") for o in options)
    # A per-option `dictionary` (input file) overriding the class default gets
    # its own column so heterogeneous option sources are unambiguous.
    include_dict = _options_need_dict_column(options, class_dictionary)

    def _append_cell(out_list: list[str], value: str, first_prefix: str, cont_prefix: str):
        text = "" if value is None else str(value)
        lines = text.splitlines()

        if not lines:
            out_list.append(first_prefix)
            return

        out_list.append(f"{first_prefix}{lines[0]}")
        for line in lines[1:]:
            out_list.append(f"{cont_prefix}{line}")

    out = []
    out.append(".. list-table::")
    if include_dict:
        # An extra Dictionary column needs widths that match the column count
        # (5 base columns + however many of path/location/dictionary are shown).
        width_map = {"Key": 14, "Path": 14, "Location": 12, "Dictionary": 16,
                     "Type": 10, "Req'd": 7, "Default": 10, "Description": 31}
        cols = ["Key"]
        if include_path:
            cols.append("Path")
        if include_location:
            cols.append("Location")
        cols.append("Dictionary")
        cols += ["Type", "Req'd", "Default", "Description"]
        out.append("   :widths: " + " ".join(str(width_map[c]) for c in cols))
    elif include_path or include_location:
        out.append("   :widths: 16 16 10 8 10 40")
    else:
        out.append("   :widths: 18 12 8 12 50")
    out.append("   :header-rows: 1")
    out.append("")
    out.append("   * - Key")
    if include_path:
        out.append("     - Path")
    if include_location:
        out.append("     - Location")
    if include_dict:
        out.append("     - Dictionary")
    out.append("     - Type")
    out.append("     - Req'd")
    out.append("     - Default")
    out.append("     - Description")

    for o in options:
        key = o.get("key", "")
        path = o.get("path", "")
        loc = o.get("location", "")
        typ = o.get("type", "")
        req = "Yes" if o.get("required", False) else "No"
        dft = "" if o.get("default", None) is None else str(o.get("default"))
        desc = o.get("description", "")

        _append_cell(out, f"``{key}``", "   * - ", "       ")
        if include_path:
            _append_cell(out, f"``{path}``" if path else "", "     - ", "       ")
        if include_location:
            _append_cell(out, f"``{loc}``" if loc else "", "     - ", "       ")
        if include_dict:
            # fall back to the class-level dictionary so every row is complete
            d = o.get("dictionary") or class_dictionary or ""
            _append_cell(out, f"``{d}``" if d else "", "     - ", "       ")
        _append_cell(out, f"``{typ}``", "     - ", "       ")
        _append_cell(out, req, "     - ", "       ")
        _append_cell(out, f"``{dft}``" if dft != "" else "", "     - ", "       ")
        _append_cell(out, desc, "     - ", "       ")

    out.append("")
    return "\n".join(out)


def _rst_usage(examples: list) -> str:
    """Render usage [{title?, comment?, snippet}] to RST."""
    if not examples:
        return "No usage example available.\n\n"
    out = []
    for ex in examples:
        title = ex.get("title")
        comment = ex.get("comment")
        snippet = (ex.get("snippet") or "").rstrip()
        if title:
            out.append(f"**{title}**")
            out.append("")
        if comment:
            out.append(comment.strip())
            out.append("")
        if snippet:
            out.append(".. code-block:: cpp")
            out.append("")
            for ln in snippet.splitlines():
                out.append("   " + ln)
            out.append("")
    return "\n".join(out)


def _opts_to_list(x) -> list[dict]:
    if not x:
        return []
    if isinstance(x, list):
        return [o for o in x if isinstance(o, dict)]
    if isinstance(x, dict):
        out = []
        for k, v in x.items():
            v = dict(v or {})
            v.setdefault("key", k)
            out.append(v)
        return out
    return []


def load_yaml_doc(yaml_path: Path) -> dict:
    """Read a *.yaml file and normalize keys."""
    raw = yaml_path.read_text(encoding="utf-8")
    if "\t" in raw:
        raw = raw.replace("\t", "    ")

    data = yaml.safe_load(raw) or {}

    usage = data.get("usage", [])
    if isinstance(usage, str):
        data["usage"] = [{"snippet": usage}] if usage.strip() else []
    elif isinstance(usage, list):
        new_usage = []
        for item in usage:
            if isinstance(item, str):
                new_usage.append({"snippet": item})
            elif isinstance(item, dict):
                new_usage.append(item)
        data["usage"] = new_usage
    else:
        data["usage"] = []

    data["options"] = _opts_to_list(data.get("options"))

    ext = data.get("externalOptions", None)
    if ext is None:
        ext = data.get("additionalOptions", None)
    data["externalOptions"] = _opts_to_list(ext)

    data.setdefault("type_name", yaml_path.stem.replace(".doc",""))
    data.setdefault("description", "")
    data.setdefault("formulation", "")
    data.setdefault("admonitions", [])
    data.setdefault("options", [])
    data.setdefault("inheritedOptions", [])
    data.setdefault("externalOptions", [])
    data.setdefault("usage", [])
    return data


def render_rst_from_yaml(y: dict, class_name: str, type_name: str | None = None) -> str:
    """
    Build an RST page from the YAML dict using canonical Sphinx API style.

    Key points:
    - First lines: document title (used by :doc: and toctree)
    - Internal headers replaced with rubrics
    - API objects declared with .. cpp:class:: (or domain of choice)
    - No hidden CSS or section headers above the title
    """

    # ---- Description block -------------------------------------------------
    description = y.get("description", "")
    description = sanitize_llm_markup(description)
    description = format_equation_for_rst(description)
    description = format_table_for_rst(description)
    description = format_code_for_rst(description)
    description = replaceInlineMath(replaceInlineReference(description))

    # ---- Formulation block -------------------------------------------------
    formulation = y.get("formulation", "")
    formulation = sanitize_llm_markup(formulation)
    formulation = format_equation_for_rst(formulation)
    formulation = format_table_for_rst(formulation)
    formulation = format_code_for_rst(formulation)
    formulation = replaceInlineMath(replaceInlineReference(formulation))

    # ---- AI provenance ------------------------------------------------------
    ai_provenance = y.get("ai_provenance") or {}

    # ---- Admonitions --------------------------------------------------------
    admonitions = y.get("admonitions") or []

    # ---- Options ------------------------------------------------------------
    inherited_options = y.get("inheritedOptions") or []
    options = y.get("options") or []
    external_options = y.get("externalOptions") or []

    # ---- Usage --------------------------------------------------------------
    usage = y.get("usage") or []

    page_title = type_name or class_name

    out = []

    # -----------------------------------------------------------------------
    # 1️⃣ Anchor (optional)
    # -----------------------------------------------------------------------
    out.append(f".. _{class_name}:\n")

    # -----------------------------------------------------------------------
    # 2️⃣ Document title (first heading in page)
    # -----------------------------------------------------------------------
    out.append(page_title)
    out.append("=" * len(page_title))
    out.append("")

    # -----------------------------------------------------------------------
    # 3️⃣ API object declaration
    # -----------------------------------------------------------------------
    # ---- Description -------------------------------------------------------
    if description.strip():
        out.extend(description.splitlines())
        out.append("")

    # ---- Formulation -------------------------------------------------------
    if formulation.strip():
        out.append(".. rubric:: Formulation")
        out.append("")
        out.extend(formulation.splitlines())
        out.append("")

    # ---- Admonitions -------------------------------------------------------
    for adm in admonitions:
        rst = _rst_admonition(adm.get("kind", "note"), adm.get("body", ""))
        out.extend(rst.splitlines())
        out.append("")

    # ---- Options dictionary (input file these options are read from) -------
    class_dictionary = y.get("dictionary") or None

    # ---- Inherited options --------------------------------------------------
    if inherited_options:
        out.append(".. rubric:: Options (inherited)")
        out.append("")
        out.extend(_rst_options_list_table(
            inherited_options, class_dictionary=class_dictionary).splitlines())
        out.append("")

    # ---- Options -----------------------------------------------------------
    if options:
        out.append(".. rubric:: Options")
        out.append("")
        # When all options share the class-level dictionary, a single note is
        # tidier than a per-row column; the table omits the column in that case.
        if class_dictionary and not _options_need_dict_column(options, class_dictionary):
            out.append(f"*Set in the* ``{class_dictionary}`` *dictionary file.*")
            out.append("")
        out.extend(_rst_options_list_table(
            options, class_dictionary=class_dictionary).splitlines())
        out.append("")

    # ---- External options --------------------------------------------------
    if external_options:
        out.append(".. rubric:: Options (external)")
        out.append("")
        out.extend(_rst_options_list_table(external_options, include_path=True).splitlines())
        out.append("")

    # ---- Usage -------------------------------------------------------------
    if usage:
        out.append(".. rubric:: Usage")
        out.append("")
        out.extend(_rst_usage(usage).splitlines())
        out.append("")

    # ---- Links -------------------------------------------------------------
    out.append(".. rubric:: Links")
    out.append("")
    out.append(f"- `Doxygen doc <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H.html>`_")
    out.append(f"- `{class_name}.H <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H_source.html>`_")
    out.append(f"- `{class_name}.C <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8C_source.html>`_")
    out.append("")

    # ---- AI provenance admonition (end of page) ----------------------------
    ai_adm = _rst_ai_provenance(ai_provenance)
    if ai_adm:
        out.extend(ai_adm.splitlines())
        out.append("")

    return "\n".join(out)




def render_rst_from_H(class_name: str, description: str, options: str, usage: str, vartable: str) -> str:
    """Build an RST page from the legacy .H-parsed sections (uses your existing transforms)."""
    # Transform description/options/usage via your existing utilities
    description = transform_content(description)
    options = transform_content(options, is_options_section=True)
    usage = transform_content(usage)
    vartable_rst = transform_vartable(vartable)

    # Title
    head = []
    head.append(f".. _{class_name}:\n")
    head.append(f"{'':=<{len(class_name)}}")
    head.append(class_name)
    head.append(f"{'':=<{len(class_name)}}\n")

    body = []
    body.append("Description")
    body.append("===========\n")
    body.append(f"{description}\n")

    if vartable_rst:
        body.append("Variables")
        body.append("=========\n")
        body.append(f"{vartable_rst}\n")

    body.append("Options")
    body.append("=======\n")
    body.append(f"{options}\n")

    body.append("Usage")
    body.append("=====\n")
    body.append(f"{usage}\n")

    body.append("Link to code")
    body.append("============\n")
    body.append(f"- `Doxygen doc <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H.html>`_")
    body.append(f"- `{class_name}.H <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H_source.html>`_")
    body.append(f"- `{class_name}.C <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8C_source.html>`_")
    body.append("")

    return "\n".join(head + body)


def generate_class_rst_files(
    ffn_lib_dirs: list[str],
    rst_output_dir: str,
    ai_content_dir: Path | None = None,
    prefer_ai: bool = False,
) -> dict:
    """
    - If a sidecar YAML exists next to a .H/.C (same stem, .yaml or .doc.yaml), use YAML.
    - (currently commented) Otherwise, parse the .H file (legacy path).
    Writes mirrored RST files into rst_output_dir and returns class_entries (ABS paths).

    ai_content_dir: root of the parallel AI-enhanced content tree (mirrors repo structure).
    prefer_ai: if True use the AI version when available, falling back to original;
               if False use the original when available, falling back to AI version.
    Both directions are silent — whichever build is running always has a complete page tree.
    """
    class_entries = {}
    out_root = Path(rst_output_dir).resolve()

    yaml_cache, mother_map = build_yaml_cache_and_mothers(ffn_lib_dirs)

    for ffn_lib_dir in ffn_lib_dirs:
        src_root = Path(ffn_lib_dir).resolve()

        # Collect canonical YAML paths from both the original and AI trees so
        # that the union covers files present in only one of the two trees.
        original_yamls: set[Path] = set()
        for p in src_root.rglob("*.yaml"):
            if not any(part in {"lnInclude", "Make"} for part in p.parts):
                original_yamls.add(p)

        ai_yamls: set[Path] = set()
        if ai_content_dir is not None:
            ai_src_root = ai_content_dir / src_root.relative_to(REPO_ROOT)
            if ai_src_root.exists():
                for p in ai_src_root.rglob("*.yaml"):
                    if not any(part in {"lnInclude", "Make"} for part in p.parts):
                        # Map back to the canonical original path for deduplication
                        orig_equiv = src_root / p.relative_to(ai_src_root)
                        ai_yamls.add(orig_equiv)

        # Process the union; resolve which file to actually read per entry.
        for yaml_path in sorted(original_yamls | ai_yamls):
            # Determine effective YAML source with symmetric fallback
            if ai_content_dir is not None:
                ai_src_root = ai_content_dir / src_root.relative_to(REPO_ROOT)
                ai_version = ai_src_root / yaml_path.relative_to(src_root)
                original_exists = yaml_path.exists()
                ai_exists = ai_version.exists()

                if prefer_ai:
                    effective = ai_version if ai_exists else yaml_path
                else:
                    effective = yaml_path if original_exists else ai_version
            else:
                effective = yaml_path

            if not effective.exists():
                continue

            class_name = yaml_path.stem.replace(".doc", "")
            spec = load_yaml_doc(effective)

            h_path = yaml_path.with_suffix(".H")
            type_name = extract_type_name_from_header(str(h_path)) if h_path.exists() else None

            cpp_class = extract_declared_class_name_from_header(str(h_path)) if h_path.exists() else class_name
            spec["inheritedOptions"] = compute_inherited_only_options(
                cpp_class,
                yaml_cache,
                mother_map
            )

            rst_text = render_rst_from_yaml(spec, class_name, type_name)

            relative_dir = yaml_path.parent.relative_to(src_root)
            output_dir = out_root.joinpath(relative_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            rst_file_path = output_dir / f"{class_name}.rst"
            write_if_changed(rst_file_path, rst_text)
            class_entries[class_name] = str(rst_file_path)

    # Remove stale RST files left over from previous runs (e.g. classes whose
    # yaml was removed, or AI-only pages when regenerating the standard tree).
    expected = {Path(p) for p in class_entries.values()}
    if out_root.exists():
        for stale in out_root.rglob("*.rst"):
            if stale not in expected:
                stale.unlink()

    return class_entries


# -----------------------------
# (Unchanged) Utility functions
# -----------------------------

def transform_vartable(content):
    content = re.sub(r'^\t|^ {4}', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\vartable', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\endvartable', '', content, flags=re.MULTILINE)
    if ("|" not in content):
        return(None)

    res  = ".. list-table:: Parameters\n"
    res += "    :widths: 50 50\n"
    res += "    :header-rows: 1\n\n"
    res += "    * - Parameter\n"
    res += "      - Description\n"
    for line in content.split('\n'):
        if ("|" in line):
            line = line.split("|")
            res += f"    * - :math:`{line[0].strip()}`\n"
            res += f"      - {line[1]}\n"

    res = replaceInlineMath(res)
    return(res)

def format_code_for_rst(transformed_content):
    return(re.sub(
        r"\\verbatim(.+?)\\endverbatim",
        lambda m: "\n.. code :: cpp\n"
                + "\n    ".join([line for line in m.group(1).splitlines()])
                + "\n",
        transformed_content,
        flags=re.DOTALL
    ))

def transform_content(content, is_options_section=False):
    """
    Kept for backward-compat with .H parsing (not used in YAML path).
    """
    content = re.sub(r'^\t|^ {4}', '', content, flags=re.MULTILINE)

    lines = content.splitlines()
    transformed_lines = []
    in_note_block = False
    in_warning_block = False
    in_info_block = False

    for i, line in enumerate(lines):
        if r"\note" in line:
            transformed_lines.append(".. note::\n\n")
            in_note_block = True
            continue
        elif r"\warning" in line:
            transformed_lines.append(".. warning::\n\n")
            in_warning_block = True
            continue
        elif is_options_section and "Parameters in" in line:
            transformed_lines.append(f".. info:: \"{line.strip()}\"\n")
            transformed_lines.append("\t<table>")
            in_info_block = True
            in_info_block_first = True
            continue

        if in_note_block:
            if line.strip() == "":
                in_note_block = False
            else:
                transformed_lines.append(f"\t{line.strip()}")
                continue

        if in_warning_block:
            if line.strip() == "":
                in_warning_block = False
            else:
                transformed_lines.append(f"\t{line.strip()}")
                continue

        if in_info_block:
            if (line.strip() == "" and in_info_block_first == False) or (i == len(lines)-1):
                transformed_lines.append(f"\t</table>\n")
                in_info_block = False
                continue

            match = re.match(r"- <b>`(.*?)`</b> - (.+)", line.strip())
            if match:
                if(in_info_block_first == False):
                    n_lines = len(transformed_lines)
                    transformed_lines[n_lines-1] += "</td>\n\t\t</tr>"
                in_info_block_first = False
                param_name = match.group(1)
                description = match.group(2).strip()
                transformed_lines.append(f"\t\t<tr>\n\t\t\t<td><strong>{param_name}</strong></td>\n\t\t\t<td>{description}")
            else:
                if(line.strip() == ""):
                    continue
                n_lines = len(transformed_lines)
                transformed_lines[n_lines-1] += (f" {line.strip()}")
            continue

        transformed_lines.append(line)

    transformed_content = "\n".join(transformed_lines)


    transformed_content = format_code_for_rst(transformed_content)
    transformed_content = format_table_for_rst(transformed_content)
    transformed_content = format_equation_for_rst(transformed_content)

    transformed_content = replaceInlineMath(transformed_content)
    transformed_content = replaceInlineReference(transformed_content)

    return transformed_content


def generate_cppapi_index(
        cppapi_folder: str,
        class_entries: dict,
        sections: list[tuple]
    ) -> None:
    """
    Parameters
    ----------
    sections: list[tuple]
        tuple = (Key word path folder in generated folder, Section name, Section Depth, Add content)
    """
    root = Path(cppapi_folder).resolve()

    # Build map: name -> path starting at 'generated/...', w/o .rst extension
    filtered_class_entries = {}
    for key, abs_path in class_entries.items():
        p = Path(abs_path).resolve()
        if "include" in p.parts:
            continue

        if "generated" in p.parts:
            idx = p.parts.index("generated")
            rel_from_generated = Path(*p.parts[idx:]).with_suffix('').as_posix()
        else:
            # Fallback: relative to cppapi_folder
            try:
                rel_from_generated = p.relative_to(root).with_suffix('').as_posix()
            except ValueError:
                rel_from_generated = p.stem  # last resort

        filtered_class_entries[key] = rel_from_generated

    cppapi_index = ".. _cppapi:\n\n"

    cppapi_index += "===============\n"
    cppapi_index += "C++ Source Code\n"
    cppapi_index += "===============\n\n"

    for sectionKey, sectionName, sectionDepth, isAddContent in sections:
        if (sectionDepth == 0):
            cppapi_index += f"{'':-<{len(sectionName)}}\n"

        cppapi_index += f"{sectionName}\n"

        if (sectionDepth == 0):
            cppapi_index += f"{'':-<{len(sectionName)}}\n\n"
        if (sectionDepth == 1):
            cppapi_index += f"{'':-<{len(sectionName)}}\n\n"
        if (sectionDepth == 2):
            cppapi_index += f"{'':~<{len(sectionName)}}\n\n"

        if (isAddContent):
            cppapi_index += ".. toctree::\n"
            cppapi_index += "   :hidden:\n"
            cppapi_index += "   :maxdepth: 1\n\n"
            for _, relpath in filtered_class_entries.items():
                if (sectionKey in relpath):
                    cppapi_index += f"   {relpath}\n"

            cppapi_index += "\n"

            cppapi_index += ".. list-table::\n"
            cppapi_index += "    :widths: 50 50\n\n"
            for name, relpath in filtered_class_entries.items():
                if (sectionKey in relpath):
                    cppapi_index += f"    * - :doc:`{relpath}`\n"

                    # Read description from ABS path
                    with open(class_entries[name], 'r', encoding='utf-8') as f:
                        for i, line in enumerate(f):
                            if (
                                i < 5
                                or line == "\n"
                                or "===" in line
                                or "Description" in line
                                or ".. _" in line
                            ):
                                continue

                            if ("Options" in line):
                                description = ""
                                break

                            description = line.replace("\n", "")
                            if ('.' in line):
                                description = line.split('.')[0]

                            break

                    cppapi_index += f"      - {description}\n"

        cppapi_index += "\n\n"

    # Add C++ Documentation YAML progress
    cppapi_index += "-----------------\n"
    cppapi_index += "C++ Documentation\n"
    cppapi_index += "-----------------\n"
    cppapi_index += "\n"
    cppapi_index += ".. toctree::\n"
    cppapi_index += "    :hidden:\n"
    cppapi_index += "    :maxdepth: 1\n"
    cppapi_index += "    yaml_progress\n"
    cppapi_index += "\n"

    cppapi_index += ":doc:`yaml_progress`\n"

    cppapi_index += "\n\n"

    write_if_changed(f"{cppapi_folder}/index.rst", cppapi_index)



#==============================================================================*
# Main function to generate both rst files and index.rst

def main():
    parser = argparse.ArgumentParser(
        description="Generate Sphinx RST files from per-class YAML docs."
    )
    parser.add_argument(
        "--ai-content-dir",
        default=None,
        metavar="DIR",
        help=(
            "Root of the parallel AI-enhanced content tree (as written by "
            "doc_ai.py --output-dir). When given, AI-enhanced YAMLs are used "
            "instead of originals (with silent fallback in both directions)."
        ),
    )
    args = parser.parse_args()

    ai_content_dir = Path(args.ai_content_dir).resolve() if args.ai_content_dir else None
    prefer_ai = ai_content_dir is not None

    # Paths to the necessary directories
    rst_output_dir = "documentation/sphinx/cppapi"  # Path to output rst files

    # Step 1: Generate rst files for all classes from YAML docs
    print("Call generate_class_rst_files (YAML) ...")
    if ai_content_dir:
        print(f"  AI content dir : {ai_content_dir}")
    os.makedirs(rst_output_dir, exist_ok=True)
    class_entries = generate_class_rst_files(
        [
            "src",
            "applications",
        ],
        rst_output_dir + "/generated",
        ai_content_dir=ai_content_dir,
        prefer_ai=prefer_ai,
    )

    # Step 2: Generate the main cppapi index
    print("Call generate_cppapi_index ...")
    generate_cppapi_index(
        cppapi_folder=rst_output_dir,
        class_entries=class_entries,
        sections=[
            # Key word, Section name, Section Depth, Add content

            # Neutronics
            ("modules/neutronics", "Neutronics", 0, True),
            ("crossSections/XS", "Cross-sections", 1, True),
            # Thermal-hydraulics
            ("modules/thermalHydraulics", "Thermal-hydraulics", 0, False),
            ("modules/thermalHydraulics/onePhase", "1-phase Solvers", 1, True),
            ("modules/thermalHydraulics/twoPhase", "2-phase Solvers", 1, True),
            ("thermophysicalProperties", "Fluid Thermophysical Properties", 1, True),
            ("porousMediaModels/phaseModels/structureModels", "Phase Models", 1, False),
            ("porousMediaModels/phaseModels/structureModels/powerModels", "Power Models", 2, True),
            ("porousMediaModels/phaseModels/structureModels/heatExchanger", "Heat Exchanger", 2, True),
            ("porousMediaModels/phaseModels/structureModels/pump", "Pump", 2, True),
            ("porousMediaModels/phaseModels/structureModels/powerOffCriterionModels", "Power Off Criterion", 2, True),
            ("porousMediaModels/phaseModels/structureModels", "Physics Models", 1, False),
            ("porousMediaModels/physicsModels/dragModels", "Drag Models", 2, True),
            ("porousMediaModels/physicsModels/heatTransferModels", "Heat transfer Models", 2, True),
            ("porousMediaModels/physicsModels/regimeMapModels", "Regime Map Models", 2, True),
            ("porousMediaModels/physicsModels/contactPartitionModels", "Contact Partition Models", 2, True),
            ("porousMediaModels/physicsModels/dispersionModels", "Dispersion Models", 2, True),
            ("porousMediaModels/physicsModels/fluidDiameterModels", "Fluid Diameter Models", 2, True),
            ("porousMediaModels/physicsModels/interfacialAreaModels", "Interfacial Area Models", 2, True),
            ("porousMediaModels/physicsModels/phaseChangeModels", "Phase Change Models", 2, True),
            ("porousMediaModels/physicsModels/phasePairs", "Phase Pairs", 2, True),
            ("porousMediaModels/physicsModels/templatedModels", "Templated Models", 2, True),
            ("porousMediaModels/physicsModels/turbulenceModels", "Turbulence Models", 2, True),
            ("porousMediaModels/physicsModels/Virtual Mass Models", "Virutal Mass Models", 2, True),
            # Offbeat lib
            ("offbeatLib", "Structural mechanics and fuel performance", 0, False),
            ("offbeatLib/physicsSubSolvers", "Physics Sub-Solvers", 1, True),
            ("offbeatLib/heatSource", "Heat Sources", 1, True),
            ("offbeatLib/fastFlux", "Fast Flux", 1, True),
            ("offbeatLib/materials/materialModel", "Materials", 1, True),
            ("offbeatLib/gapGasModel", "Gap Gas Models", 1, True),
            # Function Objects
            ("functionObjects", "Function Objects", 0, True),
            # Boundary conditions
            ("fvPatchFields", "Boundary Conditions", 0, False),
            ("fvPatchFields/neutronics", "Neutronics", 1, True),
            ("fvPatchFields/thermalHydraulics", "Thermal-hydraulics", 1, True),
            ("fvPatchFields/thermoMechanics", "Thermo-mechanics", 1, True),
            ("offbeatLib/fvPatchFields", "Thermo-mechanics", 1, True),
            # Multi-physics
            ("multiRegion", "Multi-region", 0, False),
            ("multiRegion/loop/loopModels", "Loop Models", 1, True),
            ("multiRegion/meshHandler", "Mesh Handler", 1, True),
            ("multiRegion/regionSolvers", "Region Solvers", 1, True),
            ("multiRegion/solver", "Solver", 1, True),
            # Utils
            ("profiles", "Utils", 0, False),
            ("profiles", "Profiles", 1, True),
            ("interpolationModels", "Interpolation models", 1, True),
            ("solutionControl", "Solution control", 1, True),
        ]
    )



    print("End")


#==============================================================================*
# Main

if __name__ == "__main__":
    main()

    print(f"python3 {sys.argv[0]} ... End")

#==============================================================================*