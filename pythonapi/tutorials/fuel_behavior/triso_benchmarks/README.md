# TRISO Benchmark Suite

Verification and validation cases for the OFFBEAT thermo-mechanical solver, based on
the IAEA TECDOC-1674 benchmark (Section 9). Each case simulates a 1-D spherical TRISO
fuel particle in 1-D radial geometry using the foamForNuclear Python API.

The benchmark is organised into four groups of increasing physical complexity. Each
group lives in its own subdirectory and is self-contained: it has its own geometry and
material data (`_config.py`), case builder (`case.py`), runner (`run.py`), post-processor
(`post.py`), and cleaner (`clean.py`). Automated verification tests are collected at the
root level (`test_cases_*.py`) and run with pytest.

---

## Case groups

### `cases_1_2_3` — Elastic benchmark (no irradiation)

**Physics:** purely elastic response under inner pressure and a temperature gradient.
No swelling, no creep, no time evolution.

| Case | Geometry | Loading |
|---|---|---|
| case_1 | SiC only | Thermal gradient + inner pressure |
| case_2 | IPyC only | Uniform temperature + inner pressure |
| case_3 | IPyC + SiC | Uniform temperature + inner pressure (two-layer contact) |

**Verification:** `post.py` computes the exact Lamé solution (via `analytical.py`) and
overlays it on the simulation stress profile. The tests check σ_θ at the inner face of
each layer to within 2% relative error. These cases run in seconds and carry no `slow`
marker.

---

### `cases_4` — Swelling and creep, simplified models (1000-day irradiation)

**Physics:** IPyC + SiC two-layer geometry under constant uniform temperature (1273 K)
and rising inner pressure, irradiated for 1000 days. Tests swelling and creep models
individually and in combination.

| Case | Swelling | Creep |
|---|---|---|
| case_4a | Constant isotropic rate | None |
| case_4b | None | Constant coefficient |
| case_4c | Constant isotropic | Constant coefficient |
| case_4d | Dose-dependent anisotropic (correlation A) | Constant coefficient |

**Verification:** compared to the TECDOC reference values. Cases are marked `@pytest.mark.slow`
because they run a 1000-day irradiation.

---

### `cases_5_6_7` — Full PyC model, three-layer geometry (1000-day irradiation)

**Physics:** IPyC + SiC + OPyC three-layer geometry. Realistic dose-dependent anisotropic
swelling (polynomial correlation fitted to experimental data) combined with constant creep.
No analytical solution exists; the benchmark is a code-comparison exercise.

| Case | Inner radius | Inner pressure (final) | Swelling set |
|---|---|---|---|
| case_5 | 275 µm | 15.54 MPa | A |
| case_6 | 350 µm | 26.2 MPa | A |
| case_7 | 350 µm | 26.2 MPa | B |

Cases 5 and 6 use the same swelling correlation (set A) with different fuel kernel sizes
and pressures. Case 7 uses an alternative correlation (set B) to probe sensitivity to the
PyC swelling model.

**Verification:** `post.py` checks that the final tangential stress at the inner face of
IPyC and SiC falls within the CRP-6 code-comparison band from TECDOC-1674 §9.
Tests are `@pytest.mark.slow`.

---

### `cases_8` — Temperature cycling (10 cycles, 1000 days)

**Physics:** same three-layer geometry as cases 5–7 but with a sawtooth
temperature/pressure/fluence history: temperature ramps linearly from 873 K to 1273 K
over each 100-day cycle, then drops immediately. Inner pressure and fast fluence follow
the same 10-cycle profile from TECDOC-1674. Swelling set C and a temperature-dependent
creep coefficient are used.

**Verification:** `post.py` plots σ_θ at the inner face of IPyC and SiC versus fluence,
showing the characteristic sawtooth stress evolution, and compares against digitised
TECDOC reference curves. Test is `@pytest.mark.slow`.

---

## Directory layout

```
triso_benchmarks/
├── conftest.py                # shared pytest fixtures and I/O helpers
├── test_cases_123.py          # automated tests for cases 1–3
├── test_cases_4.py            # automated tests for cases 4a–4d
├── test_cases_567.py          # automated tests for cases 5–7
├── test_cases_8.py            # automated tests for case 8
├── clean.py                   # top-level clean: delegates to groups + removes XLSXs
│
├── reference_figures/         # all TECDOC reference data (committed)
│   ├── axes.py                # axis scales for every SVG (xlim, ylim, ytick_step)
│   ├── extractSVGPlotData.py  # digitiser: SVG → XLSX
│   ├── case_5_IPyC_stress.svg
│   ├── case_5_SiC_stress.svg
│   └── …                     # one SVG per TECDOC figure; XLSXs are generated, not committed
│
├── cases_1_2_3/
│   ├── _config.py             # geometry, materials, boundary conditions
│   ├── analytical.py          # Lamé exact solution
│   ├── case.py                # builds the OpenFOAM case via the Python API
│   ├── run.py                 # runs one or more cases
│   ├── post.py                # plots simulation vs analytical
│   ├── clean.py               # removes generated case directories and figures/
│   ├── case_1/  case_2/  case_3/   # generated at runtime (not committed)
│   └── figures/               # generated at runtime (not committed)
│
├── cases_4/                   # same layout, no analytical.py
├── cases_5_6_7/               # same layout
└── cases_8/                   # same layout
```

Each `case_N/` subdirectory is an OpenFOAM case generated by `case.py`; it is not
committed to version control and is re-created by `run.py`.

---

## Reference figures

Cases 4–8 are compared against figures extracted directly from the TECDOC PDF (converted
to SVG to avoid copyright issues and to enable automated digitisation).
All reference data lives in the central `reference_figures/` directory.

### Files

| File | Status | Purpose |
|---|---|---|
| `*.svg` | committed | Primary source — the digitised TECDOC figure |
| `axes.py` | committed | Axis scales (xlim, ylim, ytick_step) for every SVG |
| `extractSVGPlotData.py` | committed | Converts SVG paths → normalised XLSX |
| `*.xlsx` | **gitignore** | Derived from SVG; auto-regenerated by `post.py` each run |

### How it works

Each `post.py` declares a `REFS` dict mapping `(case_id, layer)` to an SVG stem
(e.g. `"case_5_IPyC_stress"`). At post-processing time it:

1. Looks up the SVG in `reference_figures/<stem>.svg`.
2. Calls `extractSVGPlotData.py` to (re)generate `reference_figures/<stem>.xlsx`.
3. Reads the XLSX (`x_norm`, `y_norm` columns, normalised to [0, 1]).
4. Looks up the real axis scales from `reference_figures/axes.py` and denormalises.
5. Plots the reference curves alongside the OFFBEAT result.

### Adding a new reference figure

1. Export the TECDOC figure as SVG and place it in `reference_figures/`.
2. Add an entry to `reference_figures/axes.py` keyed by the SVG stem.
3. Add the stem to the relevant `post.py` `REFS` dict.

---

## Typical workflow

Run from inside the relevant case-group directory (all scripts use
`os.path.abspath(__file__)` so they work regardless of where Python is invoked from).

```bash
# 1. Build and run all cases in a group
cd cases_5_6_7
python run.py                  # all cases
python run.py case_5           # single case

# 2. Post-process
python post.py                 # all cases — saves figures/ inside this directory
python post.py case_6 case_7   # subset

# 3. Clean (remove generated OpenFOAM case directories)
python clean.py
```

For the elastic group (1–3), `post.py` also prints a table of σ_r and σ_θ vs the
analytical solution at each layer's inner face.

---

## Running the test suite

From the `triso_benchmarks/` root:

```bash
# Fast tests only (cases 1–3, seconds)
pytest -v test_cases_123.py

# All tests including 1000-day simulations
pytest -v

# Skip slow tests
pytest -v -m "not slow"

# Run a specific parametrised case
pytest -v test_cases_567.py -k "case_6 and IPyC"
```

The `-v` flag prints each test name and pass/fail status as it runs, which is useful
given the slow tests can take several minutes.

The tests invoke `run.py` via subprocess and cache the result within a pytest session
(running the same case twice is free). If a case has already been run manually, the test
reads the existing output directly and skips the simulation.
