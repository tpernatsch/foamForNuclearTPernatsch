"""
Analytical Lamé solutions for TRISO benchmark cases (IAEA TECDOC-1674, Ch. 9).

All three cases use pure elasticity with thermalExpansion=False, so the
analytical solution is purely mechanical (pressure loading only).

Reference values (TECDOC-1674 §9.4):
  case_1:  σ_t(r_i) =  125.19 MPa
  case_2:  σ_t(r_i) =   50.20 MPa
  case_3:  Π = −18.8 MPa  (IPyC/SiC interface pressure)
           σ_t_IPyC(r_i) =   8.8 MPa
           σ_t_SiC (r_m) = 104.4 MPa
"""

import numpy as np
from _config import CASES as _CASES, MAT as _MAT

um = 1e-6   # µm → m


# ---------------------------------------------------------------------------
# Lamé solution for a single pressurized hollow sphere
#
#   σ_r(r) = A − B/r³
#   σ_θ(r) = A + B/(2r³)
#
# where
#   A = (P_i·r_i³ − P_o·r_o³) / (r_o³ − r_i³)
#   B = (P_i − P_o)·r_i³·r_o³ / (r_o³ − r_i³)
#
# Positive stress = tensile.  Pressure BCs: σ_r(r_i)=−P_i, σ_r(r_o)=−P_o.
# ---------------------------------------------------------------------------

def _lame_AB(P_i, P_o, r_i, r_o):
    d = r_o**3 - r_i**3
    A = (P_i * r_i**3 - P_o * r_o**3) / d
    B = (P_i - P_o) * r_i**3 * r_o**3 / d
    return A, B


def sigma_r_single(r, r_i, r_o, P_i, P_o):
    A, B = _lame_AB(P_i, P_o, r_i, r_o)
    return A - B / r**3


def sigma_t_single(r, r_i, r_o, P_i, P_o):
    A, B = _lame_AB(P_i, P_o, r_i, r_o)
    return A + B / (2.0 * r**3)


# ---------------------------------------------------------------------------
# Radial displacement for a single shell (mechanical only, no thermal term)
#
#   u(r) = [(1−2ν)·A·r + (1+ν)·B/(2r²)] / E
# ---------------------------------------------------------------------------

def _u_mech(r, r_i, r_o, E, nu, P_i, P_o):
    A, B = _lame_AB(P_i, P_o, r_i, r_o)
    return ((1.0 - 2.0*nu) * A * r + (1.0 + nu) * B / (2.0 * r**2)) / E


# ---------------------------------------------------------------------------
# Interface pressure for two bonded spherical shells (analytical linear solve)
#
# Displacement continuity at r_m:
#   u_a(r_m | P_i, Π)  =  u_b(r_m | Π, P_o)
#
# Both u_a and u_b are linear in the boundary pressures, so Π can be
# found in closed form.
# ---------------------------------------------------------------------------

def interface_pressure_two_layer(
    r_ia, r_m, r_ob,
    E_a, nu_a,
    E_b, nu_b,
    P_i, P_o,
):
    """Interface pressure Π at r_m between inner shell [r_ia,r_m] and
    outer shell [r_m,r_ob], from mechanical displacement continuity.

    Returns Π [Pa]  (positive = compressive).
    """
    da = r_m**3 - r_ia**3
    db = r_ob**3 - r_m**3

    # Displacement coefficients at r_m:  u = c_A*A + c_B*B
    c_Aa = (1.0 - 2.0*nu_a) * r_m / E_a
    c_Ba = (1.0 + nu_a) / (2.0 * r_m**2 * E_a)
    c_Ab = (1.0 - 2.0*nu_b) * r_m / E_b
    c_Bb = (1.0 + nu_b) / (2.0 * r_m**2 * E_b)

    # u_a(r_m) = f_a*P_i - g_a*Π
    f_a = r_ia**3 * (c_Aa + c_Ba * r_m**3) / da
    g_a = r_m**3  * (c_Aa + c_Ba * r_ia**3) / da

    # u_b(r_m) = g_b*Π - f_b*P_o
    f_b = r_ob**3 * (c_Ab + c_Bb * r_m**3) / db
    g_b = r_m**3  * (c_Ab + c_Bb * r_ob**3) / db

    # f_a*P_i - g_a*Π = g_b*Π - f_b*P_o
    # Π = (f_a*P_i + f_b*P_o) / (g_a + g_b)
    return (f_a * P_i + f_b * P_o) / (g_a + g_b)


# ---------------------------------------------------------------------------
# Per-case analytical solution
# ---------------------------------------------------------------------------

def analytical(case_id: str, n_points: int = 200):
    """Return analytical radial and tangential stress profiles.

    Parameters
    ----------
    case_id  : "case_1", "case_2", or "case_3"
    n_points : number of radial sample points per layer

    Returns
    -------
    layers : list of dicts, one per material layer, each with keys:
        'name'    : str
        'r'       : np.ndarray [m]
        'sigma_r' : np.ndarray [Pa]
        'sigma_t' : np.ndarray [Pa]
    """
    cfg = _CASES[case_id]
    layers_cfg = cfg["layers"]
    r_i_global = cfg["inner_radius"] * um
    P_i = cfg["p_inner"]
    P_o = cfg["p_outer"]

    # Single-layer cases (case_1, case_2)
    if len(layers_cfg) == 1:
        name, r_o_um, _ = layers_cfg[0]
        r_o = r_o_um * um
        r = np.linspace(r_i_global, r_o, n_points)
        return [{
            "name": name,
            "r": r,
            "sigma_r": sigma_r_single(r, r_i_global, r_o, P_i, P_o),
            "sigma_t": sigma_t_single(r, r_i_global, r_o, P_i, P_o),
        }]

    # Two-layer case (case_3: IPyC + SiC)
    if len(layers_cfg) == 2:
        name_a, r_m_um, _ = layers_cfg[0]
        name_b, r_o_um,  _ = layers_cfg[1]
        r_m = r_m_um * um
        r_o = r_o_um * um
        mat_a = _MAT[name_a]
        mat_b = _MAT[name_b]

        Pi = interface_pressure_two_layer(
            r_ia=r_i_global, r_m=r_m, r_ob=r_o,
            E_a=mat_a["YoungModulus"], nu_a=mat_a["PoissonRatio"],
            E_b=mat_b["YoungModulus"], nu_b=mat_b["PoissonRatio"],
            P_i=P_i, P_o=P_o,
        )

        r_a = np.linspace(r_i_global, r_m, n_points)
        r_b = np.linspace(r_m, r_o, n_points)

        return [
            {
                "name": name_a,
                "r": r_a,
                "sigma_r": sigma_r_single(r_a, r_i_global, r_m, P_i, Pi),
                "sigma_t": sigma_t_single(r_a, r_i_global, r_m, P_i, Pi),
                "Pi": Pi,
            },
            {
                "name": name_b,
                "r": r_b,
                "sigma_r": sigma_r_single(r_b, r_m, r_o, Pi, P_o),
                "sigma_t": sigma_t_single(r_b, r_m, r_o, Pi, P_o),
                "Pi": Pi,
            },
        ]

    raise NotImplementedError(f"analytical() not implemented for {len(layers_cfg)}-layer case")


# ---------------------------------------------------------------------------
# Quick sanity check against TECDOC-1674 §9.4 reference values
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    for case_id, ref_str in [
        ("case_1", "ref σ_t(r_i) = 125.19 MPa"),
        ("case_2", "ref σ_t(r_i) =  50.20 MPa"),
    ]:
        res = analytical(case_id)
        r_i = _CASES[case_id]["inner_radius"] * um
        r = res[0]["r"]
        idx = np.argmin(np.abs(r - r_i))
        print(f"{case_id}: σ_t(r_i) = {res[0]['sigma_t'][idx]/1e6:8.2f} MPa  ({ref_str})")

    print("\ncase_3: IPyC + SiC")
    res = analytical("case_3")
    Pi = res[0]["Pi"]
    # Π is a contact pressure (positive=compressive); σ_r at interface = −Π
    print(f"  σ_r at interface   = {-Pi/1e6:.2f} MPa  (ref: −18.8 MPa)")
    r_i = _CASES["case_3"]["inner_radius"] * um
    r_m = _CASES["case_3"]["layers"][0][1] * um
    idx_a = np.argmin(np.abs(res[0]["r"] - r_i))
    idx_b = np.argmin(np.abs(res[1]["r"] - r_m))
    print(f"  σ_t_IPyC(r_i) = {res[0]['sigma_t'][idx_a]/1e6:7.2f} MPa  (ref:   8.8 MPa)")
    print(f"  σ_t_SiC (r_m) = {res[1]['sigma_t'][idx_b]/1e6:7.2f} MPa  (ref: 104.4 MPa)")
