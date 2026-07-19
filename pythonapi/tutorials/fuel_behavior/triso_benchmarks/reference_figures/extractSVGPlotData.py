"""
Extract all paths from an SVG file, grouped by colour.
- Resolves cumulative transforms
- Reads the legend to map hex colours to series names automatically
- Normalises coordinates against the bounding box of path id="axes"
- Outputs an Excel file (one sheet per series, columns x_norm / y_norm)
"""

import re
import sys
import math
import xml.etree.ElementTree as ET
from pathlib import Path
import pandas
import numpy as np


_NUM3  = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
_NUMS_RE = re.compile(_NUM3)

def _parse_nums(s):
    return [float(v) for v in _NUMS_RE.findall(s)]

_NAMED = {
    "red":"#ff0000","green":"#008000","blue":"#0000ff","yellow":"#ffff00",
    "cyan":"#00ffff","magenta":"#ff00ff","orange":"#ffa500","purple":"#800080",
    "white":"#ffffff","gray":"#808080","grey":"#808080","maroon":"#800000",
    "teal":"#008080","navy":"#000080","lime":"#00ff00","silver":"#c0c0c0",
    "violet":"#ee82ee","indigo":"#4b0082","pink":"#ffc0cb","brown":"#a52a2a",
}

def normalise_color(raw):
    if not raw: return None
    c = raw.strip().lower()
    if c in ("none","transparent","inherit","currentcolor"): return None
    if c in _NAMED: c = _NAMED[c]
    if re.fullmatch(r"#[0-9a-f]{3}", c): c = "#"+"".join(ch*2 for ch in c[1:])
    if re.fullmatch(r"#[0-9a-f]{6}", c): return c
    m = re.fullmatch(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", c)
    if m: return "#{:02x}{:02x}{:02x}".format(*map(int, m.groups()))
    return None

def is_black(c):  return c in {"#000000"}
def is_white(c):  return c is not None and c.lower() in {"#ffffff"}
def is_grey(c):   return c in {"#808080"}
def is_colour(c): return c is not None and not (is_white(c) or is_black(c) or is_grey(c))

def get_path_colour(attrib, style_map):
    merged = {**attrib, **style_map}
    for key in ("stroke","fill"):
        c = normalise_color(merged.get(key))
        if c and is_colour(c): return c
    return None

def has_stroke(attrib, style_map):
    """Return True if the path has a visible, non-none stroke colour."""
    merged = {**attrib, **style_map}
    c = normalise_color(merged.get("stroke", ""))
    return c is not None and is_colour(c)

def get_path_dash(attrib, style_map):
    """Return a canonical dash-array string, or '' for solid lines."""
    merged = {**attrib, **style_map}
    raw = merged.get("stroke-dasharray", "none").strip().lower()
    if raw in ("", "none"):
        return ""
    nums = _parse_nums(raw)
    return ",".join(f"{v:g}" for v in nums) if nums else ""

def parse_style(s):
    r = {}
    for d in s.split(";"):
        if ":" in d:
            k,_,v = d.partition(":")
            r[k.strip().lower()] = v.strip().lower()
    return r

def mat_id():   return [1.0,0.0,0.0,1.0,0.0,0.0]
def mat_mul(m1,m2):
    a1,b1,c1,d1,e1,f1=m1; a2,b2,c2,d2,e2,f2=m2
    return [a1*a2+c1*b2,b1*a2+d1*b2,a1*c2+c1*d2,b1*c2+d1*d2,
            a1*e2+c1*f2+e1,b1*e2+d1*f2+f1]
def mat_apply(m,x,y): a,b,c,d,e,f=m; return a*x+c*y+e, b*x+d*y+f

def parse_transform(s):
    m = mat_id()
    for fn,args_str in re.findall(r"(\w+)\s*\(([^)]*)\)",s):
        v = _parse_nums(args_str); fn = fn.lower()
        if fn=="matrix" and len(v)==6: t=v
        elif fn=="translate": t=[1,0,0,1,v[0],(v[1] if len(v)>1 else 0)]
        elif fn=="scale": sx=v[0];sy=v[1] if len(v)>1 else sx; t=[sx,0,0,sy,0,0]
        elif fn=="rotate":
            a=math.radians(v[0]); ca,sa=math.cos(a),math.sin(a)
            if len(v)==3:
                cx2,cy2=v[1],v[2]
                t=mat_mul(mat_mul([1,0,0,1,cx2,cy2],[ca,sa,-sa,ca,0,0]),[1,0,0,1,-cx2,-cy2])
            else: t=[ca,sa,-sa,ca,0,0]
        elif fn=="skewx": t=[1,0,math.tan(math.radians(v[0])),1,0,0]
        elif fn=="skewy": t=[1,math.tan(math.radians(v[0])),0,1,0,0]
        else: continue
        m = mat_mul(m,t)
    return m

def build_transform_map(root):
    transforms = {}
    def walk(elem, par):
        own = elem.attrib.get("transform","")
        mat = mat_mul(par, parse_transform(own)) if own else par
        transforms[elem] = mat
        for ch in elem: walk(ch, mat)
    walk(root, mat_id())
    return transforms

def apply_transform(pts, mat):
    return [mat_apply(mat, x, y) for x, y in pts]

def tokenize_path(d):
    for cmd_tok, num_tok in re.findall(
            r"([MmZzLlHhVvCcSsQqTtAa])|([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)", d):
        yield cmd_tok, num_tok

def extract_points(d):
    points = []; cx=cy=sx=sy=0.0; cmd,args=None,[]
    def flush():
        nonlocal cx,cy,sx,sy
        if not cmd: return
        rel=cmd.islower(); u=cmd.upper()
        def axy(dx,dy): return (cx+dx,cy+dy) if rel else (dx,dy)
        if u=="M":
            for i,(x,y) in enumerate(zip(args[0::2],args[1::2])):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
                if i==0: sx,sy=ax,ay
        elif u=="Z": cx,cy=sx,sy
        elif u=="L":
            for x,y in zip(args[0::2],args[1::2]):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
        elif u=="H":
            for x in args: ax=cx+x if rel else x; points.append((ax,cy)); cx=ax
        elif u=="V":
            for y in args: ay=cy+y if rel else y; points.append((cx,ay)); cy=ay
        elif u=="C":
            for x1,y1,x2,y2,x,y in zip(*[iter(args)]*6):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
        elif u=="S":
            for x2,y2,x,y in zip(*[iter(args)]*4):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
        elif u=="Q":
            for x1,y1,x,y in zip(*[iter(args)]*4):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
        elif u=="T":
            for x,y in zip(args[0::2],args[1::2]):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
        elif u=="A":
            for rx,ry,rot,la,sw,x,y in zip(*[iter(args)]*7):
                ax,ay=axy(x,y); points.append((ax,ay)); cx,cy=ax,ay
    for cmd_tok,num_tok in tokenize_path(d):
        if cmd_tok: flush(); cmd,args=cmd_tok,[]
        elif num_tok: args.append(float(num_tok))
    flush()
    return points

def bounding_box(pts):
    xs,ys=zip(*pts); return min(xs),min(ys),max(xs),max(ys)

def normalise(pts, x_min, y_min, x_max, y_max):
    w,h=x_max-x_min,y_max-y_min
    if w==0 or h==0: raise ValueError("Bounding box has zero width or height.")
    return [((x-x_min)/w,(y_max-y)/h) for x,y in pts]


# ---------------------------------------------------------------------------
# Legend extraction
# ---------------------------------------------------------------------------

def extract_legend(root, ns, transforms):
    """
    Returns ({hex_colour: label_string}, set_of_legend_path_ids).

    Requires a path with id="legend" in the SVG (same convention as id="axes").
    That bounding box is used to locate legend marker paths and tspan labels
    exactly.  If no id="legend" path is found, a warning is printed and empty
    results are returned.
    """

    # ── Optional: find legend bounding box (path id="legend") ───────────────
    legend_bbox = None
    for elem in root.iter(f"{ns}path"):
        if elem.attrib.get("id") == "legend":
            pts = apply_transform(
                extract_points(elem.attrib.get("d", "")), transforms[elem])
            if pts:
                legend_bbox = bounding_box(pts)
            break
    if legend_bbox:
        print(f'Legend bbox: x=[{legend_bbox[0]:.1f}, {legend_bbox[2]:.1f}]'
              f'  y=[{legend_bbox[1]:.1f}, {legend_bbox[3]:.1f}]')
    else:
        print('WARNING: no path with id="legend" found — legend cannot be extracted. '
              'Add a path with id="legend" to the SVG to enable legend detection.')
        return {}, set()

    # ── Collect all tspan text in global SVG coordinates ────────────────────
    # tspan x/y are LOCAL; apply accumulated transforms to compare against
    # marker centroids (which are already in global coordinates).
    all_labels = []  # (global_y, global_x, text)
    for elem in root.iter(f"{ns}tspan"):
        text = (elem.text or "").strip()
        if not text: continue
        if re.fullmatch(r"[+-]?\d+\.?\d*", text): continue  # skip tick numbers
        y_str = elem.attrib.get("y")
        x_str = elem.attrib.get("x")
        if y_str is None: continue
        try:
            lx = float(x_str) if x_str else 0.0
            ly = float(y_str)
            gx, gy = mat_apply(transforms.get(elem, mat_id()), lx, ly)
            all_labels.append((gy, gx, text))
        except ValueError:
            pass

    if not all_labels:
        return {}, set()

    # ── Build label_ys: {global_y: text} for legend labels only ─────────────
    lb_xmin, lb_ymin, lb_xmax, lb_ymax = legend_bbox
    label_ys = {gy: text for gy, gx, text in all_labels
                if lb_xmin <= gx <= lb_xmax and lb_ymin <= gy <= lb_ymax}

    if not label_ys:
        return {}, set()

    # ── Find legend marker paths and populate legend_path_ids ────────────────
    legend_path_ids = set()
    best_marker = {}  # colour -> (npts, mean_y, pid)

    lb_xmin, lb_ymin, lb_xmax, lb_ymax = legend_bbox
    for i, elem in enumerate(root.iter(f"{ns}path")):
        if elem.attrib.get("id") == "legend":
            continue  # skip the legend bbox rectangle itself
        attrib    = {k.lower(): v.lower() for k, v in elem.attrib.items()}
        style_map = parse_style(attrib.get("style",""))
        colour    = get_path_colour(attrib, style_map)
        if not colour: continue

        d     = elem.attrib.get("d","")
        t_pts = apply_transform(extract_points(d), transforms[elem])
        if not t_pts: continue
        npts  = len(t_pts)

        x_cen = sum(p[0] for p in t_pts) / npts
        y_cen = sum(p[1] for p in t_pts) / npts

        if not (lb_xmin <= x_cen <= lb_xmax and lb_ymin <= y_cen <= lb_ymax):
            continue  # outside legend bbox — not a legend element
        pid = elem.attrib.get("id", f"path_{i}")
        legend_path_ids.add(pid)          # exclude from data regardless of npts
        # Keep the representative marker for this colour: prefer fewer points,
        # then higher up in the legend (smaller y = earlier row).
        if colour not in best_marker or (npts, y_cen) < best_marker[colour][:2]:
            best_marker[colour] = (npts, y_cen, pid)

    if not best_marker:
        return {}, legend_path_ids

    # ── Sort both by global y and zip to assign names by rank ────────────────
    sorted_labels  = sorted(label_ys.items())
    sorted_markers = sorted(best_marker.items(), key=lambda kv: kv[1][1])

    if len(sorted_labels) != len(sorted_markers):
        n = min(len(sorted_labels), len(sorted_markers))
        print(f"  Warning: {len(sorted_markers)} marker(s) but {len(sorted_labels)} label(s) — "
              f"matching first {n}")
        sorted_labels  = sorted_labels[:n]
        sorted_markers = sorted_markers[:n]

    colour_to_label = {}  # colour -> name
    for (_, name), (colour, _) in zip(sorted_labels, sorted_markers):
        colour_to_label[colour] = name

    return colour_to_label, legend_path_ids


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------

def extract_all_colours(svg_path):
    tree = ET.parse(svg_path)
    root = tree.getroot()
    ns_match = re.match(r"\{.*?\}", root.tag)
    ns = ns_match.group(0) if ns_match else ""
    transforms = build_transform_map(root)

    # Axes bounding box
    axes_bbox = None
    for elem in root.iter(f"{ns}path"):
        if elem.attrib.get("id") == "axes":
            pts = apply_transform(extract_points(elem.attrib.get("d","")), transforms[elem])
            if not pts: raise ValueError('Path id="axes" has no points.')
            axes_bbox = bounding_box(pts)
            break
    if axes_bbox is None:
        raise ValueError('No path with id="axes" found.')
    x_min,y_min,x_max,y_max = axes_bbox
    print(f'Axes bbox: x=[{x_min:.4f}, {x_max:.4f}]  y=[{y_min:.4f}, {y_max:.4f}]')

    # Legend colour→name mapping + set of path IDs to skip
    colour_names, legend_path_ids = extract_legend(root, ns, transforms)
    if colour_names:
        print("Legend detected:")
        for c, n in colour_names.items():
            print(f"  {c}  ->  '{n}'")
    else:
        print("No legend found — using hex codes as series names.")
    print(f"  Skipping {len(legend_path_ids)} legend path(s)")

    # Collect data paths, keyed by colour only.
    # Dash pattern is intentionally NOT used as a key here: some SVG exporters
    # apply stroke-dasharray inconsistently across segments of the same logical
    # curve, which would split one curve into incomplete sub-series and produce
    # zigzag artefacts.  Dash information is only used in extract_legend (for
    # legend label matching); the data itself is always grouped by colour.
    # Two-pass collection: prefer stroked paths.  Some SVG renderers represent
    # dashed lines as sequences of tiny filled rectangles (stroke=none).  These
    # must be captured when no stroked paths exist for that colour; but when a
    # colour has BOTH stroked data segments AND filled scatter markers, collecting
    # only the stroked segments avoids zigzag artefacts from mixing the two.
    colour_paths      = {}   # colour -> {pid: [pts]}  (stroked only)
    colour_paths_fill = {}   # colour -> {pid: [pts]}  (fill-only, fallback)

    for i, elem in enumerate(root.iter(f"{ns}path")):
        pid = elem.attrib.get("id", f"path_{i}")
        if pid in legend_path_ids:
            continue

        attrib    = {k.lower(): v.lower() for k, v in elem.attrib.items()}
        style_map = parse_style(attrib.get("style",""))
        colour    = get_path_colour(attrib, style_map)
        if not colour: continue

        stroked = has_stroke(attrib, style_map)

        raw_pts   = extract_points(elem.attrib.get("d",""))
        world_pts = apply_transform(raw_pts, transforms[elem])
        if not world_pts: continue

        _eps = (x_max - x_min) * 0.01
        _cx  = sum(p[0] for p in world_pts) / len(world_pts)
        _cy  = sum(p[1] for p in world_pts) / len(world_pts)
        if not ((x_min - _eps) <= _cx <= (x_max + _eps)
                and (y_min - _eps) <= _cy <= (y_max + _eps)):
            continue

        _TOL = 0.02
        norm_pts = [(x, y) for x, y in normalise(world_pts, *axes_bbox)
                    if -_TOL <= x <= 1 + _TOL and -_TOL <= y <= 1 + _TOL]
        if not norm_pts:
            continue

        # Skip grid lines: axis-aligned paths spanning most of the plot area.
        _xs = [p[0] for p in norm_pts]
        _ys = [p[1] for p in norm_pts]
        _xspan = max(_xs) - min(_xs)
        _yspan = max(_ys) - min(_ys)
        if (_yspan < 0.015 and _xspan > 0.50) or (_xspan < 0.015 and _yspan > 0.50):
            continue

        if stroked:
            colour_paths.setdefault(colour, {}).setdefault(pid, []).extend(norm_pts)
        else:
            colour_paths_fill.setdefault(colour, {}).setdefault(pid, []).extend(norm_pts)

    # Merge fill-only colours that have no stroked counterpart (dashed-line
    # curves rendered as filled rectangles, e.g. Japan in some TECDOC figures).
    for colour, paths in colour_paths_fill.items():
        if colour not in colour_paths:
            colour_paths[colour] = paths

    if not colour_paths:
        print("No coloured data paths found.")
        return {}, {}

    return colour_paths, colour_names


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if not (2 <= len(sys.argv) <= 3):
        print("Usage: python extractSVGPlotData.py <file.svg> [out_dir]")
        sys.exit(1)
    svg_file = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) == 3 else "."

    colour_paths, colour_names = extract_all_colours(svg_file)
    print(f'Found {len(colour_paths)} colour(s)')

    # Merge all segments of each colour and sort by x_norm so that the curve
    # is reconstructed correctly regardless of SVG path element order.
    # Sorting also makes each series monotone in x, which is correct for
    # fluence-vs-stress benchmark plots.
    def _label(c):
        return colour_names.get(c, c.lstrip("#"))

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    excel_out = out_path / f"{Path(svg_file).stem}.xlsx"
    total = 0
    with pandas.ExcelWriter(excel_out) as f:
        for colour, paths in colour_paths.items():
            all_pts = sorted(
                [pt for pts_list in paths.values() for pt in pts_list],
                key=lambda p: p[0]
            )
            df = pandas.DataFrame(all_pts, columns=["x_norm", "y_norm"])
            total += len(df)
            sheet = re.sub(r'[\\/*?\[\]:]', '_', _label(colour))[:31]
            df.to_excel(f, sheet_name=sheet, index=False)
    print(f"Wrote {total} points across {len(colour_paths)} series → {excel_out}")