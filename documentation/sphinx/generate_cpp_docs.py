"""
Generate the doc for Sphinx RST files

Now reads per-class YAML docs placed next to the .H/.C files:
  e.g.  myClass.H, myClass.C, myClass.yaml

YAML fields (all optional except type_name is recommended):
  type_name: str
  summary: str (can be multiline with | )
  admonitions:
    - { kind: warning|note|info|tip|..., body: str }
  options:
    - { key: str, type: str, required: bool, default: any, description: str }
  usage:
    - { title?: str, comment?: str, snippet: str }

The rest of the script (release notes, index, etc.) is unchanged.
"""

#==============================================================================*
# Imports

import os
import re
import sys
from pathlib import Path

# NEW: YAML
import yaml

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
        text = text.replace(f"${refValue}$", f":math:`{refValue}`")

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
    tableRow = []
    textSplit = text.split('\n')
    isRead = False
    new_text = text
    for i, line in enumerate(textSplit):
        if ((r"\f[" in line or r"$$" in line) and not isRead):
            textSplit[i] = ".. math::"
            isRead = True
            continue
        if ((r"\f]" in line or r"$$" in line) and isRead):
            isRead = False

            for idx in tableRow:
                textSplit[idx] = f"    {textSplit[idx].strip()}"

            textSplit[i] = ""
            new_text = '\n'.join(textSplit)

            tableRow = []
            continue

        if (isRead):
            tableRow.append(i)

    return(new_text)

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


# ============================================================
# NEW: YAML-based generation (replaces .H parsing in practice)
# ============================================================

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


def _rst_options_list_table(options: list) -> str:
    """Render options (key/type/required/default/description) as a list-table."""
    if not options:
        return "No options list available.\n\n"
    out = []
    out.append(".. list-table::")
    out.append("   :widths: 18 12 8 12 50")
    out.append("   :header-rows: 1")
    out.append("")
    out.append("   * - Key")
    out.append("     - Type")
    out.append("     - Req'd")
    out.append("     - Default")
    out.append("     - Description")
    for o in options:
        key = o.get("key","")
        typ = o.get("type","")
        req = "Yes" if o.get("required", False) else "No"
        dft = "" if o.get("default", None) is None else str(o.get("default"))
        desc = o.get("description","")
        out.append(f"   * - ``{key}``")
        out.append(f"     - ``{typ}``")
        out.append(f"     - {req}")
        out.append(f"     - ``{dft}``" if dft != "" else "     - ")
        out.append(f"     - {desc}")
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


def load_yaml_doc(yaml_path: Path) -> dict:
    """Read a *.yaml file and normalize keys."""
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    data.setdefault("type_name", yaml_path.stem.replace(".doc",""))
    data.setdefault("description", "")
    data.setdefault("admonitions", [])
    data.setdefault("options", [])
    data.setdefault("usage", [])
    return data


# def render_rst_from_yaml(y: dict, class_name: str) -> str:
#     """Build an RST page from the YAML dict using Sphinx object markup."""

#     # ---- Description block -------------------------------------------------
#     summary = y.get("description", "")
#     summary = format_equation_for_rst(summary)
#     summary = format_table_for_rst(summary)
#     summary = format_code_for_rst(summary)
#     summary = replaceInlineMath(replaceInlineReference(summary))

#     # ---- Admonitions --------------------------------------------------------
#     admonitions = y.get("admonitions") or []

#     # ---- Options ------------------------------------------------------------
#     options = y.get("options") or []

#     # ---- Usage --------------------------------------------------------------
#     usage = y.get("usage") or []

#     # -----------------------------------------------------------------------
#     # Anchor (safe to keep; does not affect TOC)
#     # -----------------------------------------------------------------------
#     out = []
#     out.append(f".. _{class_name}:\n")

#     # -----------------------------------------------------------------------
#     # Class declaration (THIS replaces the ===== title)
#     # -----------------------------------------------------------------------
#     out.append(f".. cpp:class:: {class_name}\n")

#     # Everything that follows must be indented to belong to the class
#     indent = "   "

#     # ---- Description -------------------------------------------------------
#     if summary.strip():
#         for line in summary.splitlines():
#             out.append(indent + line)
#         out.append("")

#     # ---- Admonitions --------------------------------------------------------
#     for adm in admonitions:
#         rst = _rst_admonition(adm.get("kind", "note"), adm.get("body", ""))
#         for line in rst.splitlines():
#             out.append(indent + line)
#         out.append("")

#     # ---- Options ------------------------------------------------------------
#     if options:
#         out.append(indent + ".. rubric:: Options\n")
#         table = _rst_options_list_table(options)
#         for line in table.splitlines():
#             out.append(indent + line)
#         out.append("")

#     # ---- Usage --------------------------------------------------------------
#     if usage:
#         out.append(indent + ".. rubric:: Usage\n")
#         usage_rst = _rst_usage(usage)
#         for line in usage_rst.splitlines():
#             out.append(indent + line)
#         out.append("")

#     # ---- Links --------------------------------------------------------------
#     out.append(indent + ".. rubric:: Links\n")
#     out.append(indent + f"- `Doxygen doc <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H.html>`_")
#     out.append(indent + f"- `{class_name}.H <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H_source.html>`_")
#     out.append(indent + f"- `{class_name}.C <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8C_source.html>`_")
#     out.append("")

#     return "\n".join(out)

def render_rst_from_yaml(y: dict, class_name: str) -> str:
    """
    Build an RST page from the YAML dict using canonical Sphinx API style.

    Key points:
    - First lines: document title (used by :doc: and toctree)
    - Internal headers replaced with rubrics
    - API objects declared with .. cpp:class:: (or domain of choice)
    - No hidden CSS or section headers above the title
    """

    # ---- Description block -------------------------------------------------
    summary = y.get("description", "")
    summary = format_equation_for_rst(summary)
    summary = format_table_for_rst(summary)
    summary = format_code_for_rst(summary)
    summary = replaceInlineMath(replaceInlineReference(summary))

    # ---- Admonitions --------------------------------------------------------
    admonitions = y.get("admonitions") or []

    # ---- Options ------------------------------------------------------------
    options = y.get("options") or []

    # ---- Usage --------------------------------------------------------------
    usage = y.get("usage") or []

    out = []

    # -----------------------------------------------------------------------
    # 1️⃣ Anchor (optional)
    # -----------------------------------------------------------------------
    out.append(f".. _{class_name}:\n")

    # -----------------------------------------------------------------------
    # 2️⃣ Document title (first heading in page)
    # -----------------------------------------------------------------------
    out.append(class_name)
    out.append("=" * len(class_name))
    out.append("")

    # -----------------------------------------------------------------------
    # 3️⃣ API object declaration
    # -----------------------------------------------------------------------
    out.append(f".. cpp:class:: {class_name}\n")
    indent = "   "

    # ---- Description -------------------------------------------------------
    if summary.strip():
        for line in summary.splitlines():
            out.append(indent + line)
        out.append("")

    # ---- Admonitions -------------------------------------------------------
    for adm in admonitions:
        rst = _rst_admonition(adm.get("kind", "note"), adm.get("body", ""))
        for line in rst.splitlines():
            out.append(indent + line)
        out.append("")

    # ---- Options -----------------------------------------------------------
    if options:
        out.append(indent + ".. rubric:: Options\n")
        table = _rst_options_list_table(options)
        for line in table.splitlines():
            out.append(indent + line)
        out.append("")

    # ---- Usage -------------------------------------------------------------
    if usage:
        out.append(indent + ".. rubric:: Usage\n")
        usage_rst = _rst_usage(usage)
        for line in usage_rst.splitlines():
            out.append(indent + line)
        out.append("")

    # ---- Links -------------------------------------------------------------
    out.append(indent + ".. rubric:: Links\n")
    out.append(indent + f"- `Doxygen doc <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H.html>`_")
    out.append(indent + f"- `{class_name}.H <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H_source.html>`_")
    out.append(indent + f"- `{class_name}.C <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8C_source.html>`_")
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


def generate_class_rst_files(ffn_lib_dirs: list[str], rst_output_dir: str) -> dict:
    """
    - If a sidecar YAML exists next to a .H/.C (same stem, .yaml or .doc.yaml), use YAML.
    - (currently commented) Otherwise, parse the .H file (legacy path).
    Writes mirrored RST files into rst_output_dir and returns class_entries (ABS paths).
    """
    class_entries = {}
    out_root = Path(rst_output_dir).resolve()

    for ffn_lib_dir in ffn_lib_dirs:
        src_root = Path(ffn_lib_dir).resolve()

        # 1) Prefer YAML where available
        for yaml_path in src_root.rglob("*.yaml"):
            if any(p in {"lnInclude", "Make"} for p in yaml_path.parts):
                continue

            class_name = yaml_path.stem.replace(".doc", "")
            spec = load_yaml_doc(yaml_path)
            rst_text = render_rst_from_yaml(spec, class_name)

            relative_dir = yaml_path.parent.relative_to(src_root)
            output_dir = out_root.joinpath(relative_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            rst_file_path = output_dir / f"{class_name}.rst"
            rst_file_path.write_text(rst_text, encoding="utf-8")
            class_entries[class_name] = str(rst_file_path)

        # # 2) For any .H without a YAML sibling, fall back to legacy parsing
        # for h_path in src_root.rglob("*.H"):
        #     if any(p in {"lnInclude", "Make"} for p in h_path.parts):
        #         continue

        #     # Skip if RST already generated from YAML for this class
        #     class_name = h_path.stem
        #     if class_name in class_entries:
        #         continue

        #     # Look for sidecar YAML (.yaml or .doc.yaml); if found, it would have been handled above
        #     yaml_sidecars = [
        #         h_path.with_suffix(".yaml"),
        #         h_path.with_suffix(".doc.yaml"),
        #     ]
        #     if any(y.exists() for y in yaml_sidecars):
        #         continue

        #     # Legacy .H parsing
        #     description, options, usage, vartable = extract_descriptions_and_usage(str(h_path))
        #     rst_text = render_rst_from_H(class_name, description, options, usage, vartable)

        #     relative_dir = h_path.parent.relative_to(src_root)
        #     output_dir = out_root.joinpath(relative_dir)
        #     output_dir.mkdir(parents=True, exist_ok=True)

        #     rst_file_path = output_dir / f"{class_name}.rst"
        #     rst_file_path.write_text(rst_text, encoding="utf-8")
        #     class_entries[class_name] = str(rst_file_path)

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

    with open(f"{cppapi_folder}/index.rst", 'w', encoding='utf-8') as f:
        f.write(cppapi_index)



from pathlib import Path
import os

def append_toctree_for_folder_direct_subfolders(target_rst: str, folder: str, maxdepth: int = 1) -> None:
    """
    Append a Sphinx toctree block listing .rst files that are exactly one level
    below `folder` (i.e., direct subfolders only; no nesting).
    """
    folder_path = Path(folder).resolve()

    # Collect *.rst files in *direct* subfolders only (depth = 1)
    rst_files = []
    for f in folder_path.rglob("*.rst"):
        try:
            parts = f.relative_to(folder_path).parts
        except ValueError:
            # Shouldn't happen because f is under folder_path, but just in case
            continue
        # We want exactly 2 parts: ('subfolder', 'file.rst')
        if len(parts) == 2 and f.is_file():
            rst_files.append(f)

    rst_files = sorted(rst_files)

    if not rst_files:
        print(f"No .rst files found in direct subfolders of {folder_path}")
        return

    # Build toctree block
    toctree_block = [
        "\n.. toctree::",
        f"   :maxdepth: {maxdepth}",
        "",
    ]
    for f in rst_files:
        rel_path = os.path.relpath(f, start=Path(target_rst).parent)
        toctree_block.append(f"   {rel_path}")

    with open(target_rst, "a", encoding="utf-8") as out_file:
        out_file.write("\n".join(toctree_block) + "\n")

    print(f"Appended toctree (direct subfolders only) with {len(rst_files)} entries to {target_rst}")


def insert_toctree_at_placeholder_for_folder_direct_subfolders(
    target_rst: str,
    folder: str,
    maxdepth: int = 1,
    placeholder: str = "toctreeHere",
    append_if_missing: bool = True,
) -> None:
    """
    Replace a line containing `placeholder` with a Sphinx toctree listing *.rst files
    that are exactly one level below `folder` (direct subfolders only).

    - Preserves the placeholder line's leading indentation for the directive.
    - Ensures a blank line before the directive (adds one if the previous line isn't blank).
    - Indents options and entries by exactly 3 spaces relative to the directive.
    - Uses forward slashes in entries for Sphinx portability.
    """
    target_path = Path(target_rst).resolve()
    folder_path = Path(folder).resolve()

    # Collect *.rst files in direct subfolders only (depth = 1)
    rst_files = []
    for f in folder_path.rglob("*.rst"):
        try:
            parts = f.relative_to(folder_path).parts
        except ValueError:
            continue
        if len(parts) == 2 and f.is_file():
            rst_files.append(f)

    rst_files = sorted(rst_files)
    if not rst_files:
        print(f"[info] No .rst files found in direct subfolders of {folder_path}")
        return

    # Read target file
    try:
        with open(target_path, "r", encoding="utf-8") as fh:
            content = fh.read()
    except FileNotFoundError:
        print(f"[error] Target RST not found: {target_path}")
        return

    # Work line-by-line to precisely control indentation and blank lines
    lines = content.splitlines(keepends=False)

    # Find placeholder line index and its indent
    placeholder_idx = None
    placeholder_indent = ""
    for i, line in enumerate(lines):
        # Match a line that contains only optional whitespace + placeholder + optional whitespace
        if re.fullmatch(rf"\s*{re.escape(placeholder)}\s*", line):
            placeholder_idx = i
            placeholder_indent = re.match(r"\s*", line).group(0)
            break

    # Build toctree with correct indentation rules
    INDENT = "   "  # 3 spaces per Sphinx convention
    directive_indent = placeholder_indent  # directive starts at the placeholder's indent
    opt_indent = directive_indent + INDENT
    entry_indent = directive_indent + INDENT

    toctree_lines = [
        f"{directive_indent}.. toctree::",
        f"{opt_indent}:maxdepth: {maxdepth}",
        f"{directive_indent}",  # blank line separating options from entries (same indent as directive)
    ]

    for f in rst_files:
        rel_path = os.path.relpath(f, start=target_path.parent).replace(os.sep, "/")
        # You can drop the .rst suffix if you prefer; Sphinx accepts both.
        toctree_lines.append(f"{entry_indent}{rel_path}")

    # Ensure there is a blank line BEFORE the directive
    # If the placeholder is at index k, ensure lines[k-1] is blank (or it's the top of file)
    def ensure_blank_line_before(lines_list, insert_at_index, indent_for_blank=""):
        if insert_at_index <= 0:
            # Insert a blank line at the very top only if the first line isn't already blank
            return lines_list
        prev = lines_list[insert_at_index - 1]
        if prev.strip() != "":
            # Insert a truly blank line (no spaces) as Sphinx blank line
            lines_list.insert(insert_at_index, "")
        return lines_list

    if placeholder_idx is not None:
        # Replace the placeholder line with toctree_lines (injecting a blank line before if needed)
        lines = ensure_blank_line_before(lines, placeholder_idx)
        # After possibly inserting, the placeholder_idx may shift by +1
        # Recompute if a blank line was inserted
        # We detect by checking if the line at placeholder_idx is empty now
        if lines[placeholder_idx].strip() == "":
            placeholder_idx += 1

        # Remove the placeholder line
        del lines[placeholder_idx]

        # Insert toctree block at the placeholder position
        for j, tline in enumerate(toctree_lines):
            lines.insert(placeholder_idx + j, tline)

        new_content = "\n".join(lines) + "\n"
        with open(target_path, "w", encoding="utf-8") as fh:
            fh.write(new_content)

        print(
            f"[ok] Inserted toctree (direct subfolders only, {len(rst_files)} entries) "
            f"at placeholder '{placeholder}' in {target_path}"
        )
        return

    # If placeholder not found
    warn_msg = f"[warn] Placeholder '{placeholder}' not found in {target_path}."
    if append_if_missing:
        # Append with a preceding blank line for safety
        if lines and lines[-1].strip() != "":
            lines.append("")
        # For appending, left-align the directive at column 0
        directive_indent = ""
        opt_indent = "   "
        entry_indent = "   "
        toctree_lines = [
            f"{directive_indent}.. toctree::",
            f"{opt_indent}:maxdepth: {maxdepth}",
            f"{directive_indent}",
        ]
        for f in rst_files:
            rel_path = os.path.relpath(f, start=target_path.parent).replace(os.sep, "/")
            toctree_lines.append(f"{entry_indent}{rel_path}")

        lines.extend(toctree_lines)
        new_content = "\n".join(lines) + "\n"
        with open(target_path, "w", encoding="utf-8") as fh:
            fh.write(new_content)
        print(warn_msg + " Appended toctree at end of file instead.")
    else:
        print(warn_msg + " No changes made.")




def append_toctree_for_folder_recursive(target_rst: str, folder: str, maxdepth: int = 1) -> None:
    """
    Append a Sphinx toctree block at the end of `target_rst` file,
    listing all .rst files found in subfolders of `folder` (excluding root-level files).

    Parameters
    ----------
    target_rst : str
        Path to the RST file where the toctree will be appended.
    folder : str
        Path to the folder containing .rst files (e.g., dragModels).
    maxdepth : int
        Depth for the toctree (default = 1).
    """
    folder_path = Path(folder).resolve()
    # Only include .rst files that are NOT directly under folder_path
    rst_files = sorted([f for f in folder_path.rglob("*.rst") if f.parent != folder_path])

    if not rst_files:
        print(f"No .rst files found in subfolders of {folder_path}")
        return

    # Build toctree block
    toctree_block = [
        "\n.. toctree::",
        f"   :maxdepth: {maxdepth}",
        "",
    ]
    for f in rst_files:
        rel_path = os.path.relpath(f, start=Path(target_rst).parent)
        toctree_block.append(f"   {rel_path}")

    # Append to target RST file
    with open(target_rst, "a", encoding="utf-8") as out_file:
        out_file.write("\n".join(toctree_block) + "\n")

    print(f"Appended recursive toctree with {len(rst_files)} entries to {target_rst}")





#==============================================================================*
# Main function to generate both rst files and index.rst

def main():
    # Paths to the necessary directories
    rst_output_dir = "documentation/sphinx/cppapi"  # Path to output rst files

    # Step 1: Generate rst files for all classes from YAML docs
    print("Call generate_class_rst_files (YAML) ...")
    os.makedirs(rst_output_dir, exist_ok=True)
    class_entries = generate_class_rst_files(
        [
            "src",
            "applications",
        ],
        rst_output_dir+"/generated"
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

    # Step 3: Append indexes in dynamic files
    insert_toctree_at_placeholder_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/GeN-Foam/coupling.rst",
        folder="documentation/sphinx/cppapi/generated/multiRegion/loop/loopModels",
        maxdepth=1
    )
    insert_toctree_at_placeholder_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/initialAndBC.rst",
        folder="documentation/sphinx/cppapi/generated/fvPatchFields/thermalHydraulics/",
        maxdepth=1
    )
    insert_toctree_at_placeholder_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/neutronics/initialAndBC.rst",
        folder="documentation/sphinx/cppapi/generated/fvPatchFields/neutronics/",
        maxdepth=1
    )
    insert_toctree_at_placeholder_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/cppapi/generated/porousMediaModels/phaseModels/structureModels/powerModels/powerModel.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/phaseModels/structureModels/powerModels/",
        maxdepth=1
    )
    insert_toctree_at_placeholder_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/cppapi/generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/powerOffCriterionModel.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/turbulenceProperties.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/turbulenceModels/",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/FFdrag.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/FSdrag.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/multipliersDrag.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/contactPartitionModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/contactPartitionModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/dispersionModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/dispersionModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/twoPhaseSpecific/fluidDiameterModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/fluidDiameterModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/FFHeatTransferCoefficientModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/FSHeatTransferCoefficientModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/interfacialAreaModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/interfacialAreaModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/regimeMapModels/index.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/regimeMapModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/virtualMassCoefficientModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/virtualMassModels/",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/phaseChangeModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/phaseChangeModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/latentHeatModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/latentHeatModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/phaseProperties/physicsModels/saturationModels.rst",
        folder="documentation/sphinx/cppapi/generated/porousMediaModels/physicsModels/saturationModels",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/thermalHydraulics/thermalHydraulicsSolvers.rst",
        folder="documentation/sphinx/cppapi/generated/modules/thermalHydraulics/",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/modules/neutronics/neutronicsSolvers.rst",
        folder="documentation/sphinx/cppapi/generated/modules/neutronics/",
        maxdepth=1
    )
    append_toctree_for_folder_direct_subfolders(
        target_rst="documentation/sphinx/usersguide/GeN-Foam/postProcessing.rst",
        folder="documentation/sphinx/cppapi/generated/functionObjects/",
        maxdepth=1
    )

    print("End")


#==============================================================================*
# Main

if __name__ == "__main__":
    main()

    print(f"python3 {sys.argv[0]} ... End")

#==============================================================================*
