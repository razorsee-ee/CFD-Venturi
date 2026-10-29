"""Hypothetical pyrolysis gas for the CFD case: composition, mixture properties, operating point.

Flows are t/h. Every hydrocarbon cut is represented by one surrogate molecule.
Condensation is ignored: at 470 C everything is treated as vapour.
"""
import math

R = 8314.46  # J/(kmol K)

# name: (t/h, surrogate, M kg/kmol, cp at ~400-470 C J/(kg K), mu Pa s, k W/(m K))
COMPONENTS = {
    "N2":           (12 * 0.20, "N2",      28.0, 1120, 3.4e-5, 0.052),
    "CO2":          (12 * 0.10, "CO2",     44.0, 1150, 3.2e-5, 0.050),
    "light HC":     (12 * 0.70, "C2H6",    30.1, 3200, 1.9e-5, 0.075),
    "water vapour": (9,         "H2O",     18.0, 2100, 2.6e-5, 0.055),
    "gasoline":     (6,         "C7H16",  100.2, 3100, 1.4e-5, 0.045),
    "mazut":        (33,        "C18H38", 254.5, 3000, 1.0e-5, 0.035),
    "heavy oil":    (16,        "C25H52", 352.7, 3000, 0.9e-5, 0.030),
}

T_GAS = 470 + 273.15    # K
P = 101325.0            # Pa, assumed ~atmospheric at the scrubber outlet
D_IN, D_THROAT = 0.884, 0.560

OIL_TPH = 32.0          # spray: heavy oil
T_OIL = 100 + 273.15
OIL = dict(rho=1150.0, cp=1950.0, mu=0.046, sigma=0.030)  # at 100 C: 1.2 kg/l cold, 40 cSt
D_ORIFICE = 0.020       # nozzle orifice, assumed


def main():
    m_tot = sum(c[0] for c in COMPONENTS.values())
    n = {k: c[0] * 1000 / c[2] for k, c in COMPONENTS.items()}  # kmol/h
    n_tot = sum(n.values())
    M = m_tot * 1000 / n_tot
    w = {k: c[0] / m_tot for k, c in COMPONENTS.items()}
    x = {k: n[k] / n_tot for k in n}
    cp = sum(w[k] * c[3] for k, c in COMPONENTS.items())
    # Wilke-like shortcut: mole-fraction weighting with sqrt(M) (Herning-Zipperer)
    sm = {k: math.sqrt(c[2]) for k, c in COMPONENTS.items()}
    mu = sum(x[k] * c[4] * sm[k] for k, c in COMPONENTS.items()) / sum(x[k] * sm[k] for k in x)
    kappa = sum(x[k] * c[5] * sm[k] for k, c in COMPONENTS.items()) / sum(x[k] * sm[k] for k in x)
    gamma = cp / (cp - R / M)

    print(f"{'component':14s} {'t/h':>6s} {'surrogate':>9s} {'M':>6s} {'w,mass':>7s} {'x,mol':>7s}")
    for k, c in COMPONENTS.items():
        print(f"{k:14s} {c[0]:6.1f} {c[1]:>9s} {c[2]:6.1f} {w[k]:7.3f} {x[k]:7.3f}")
    print(f"{'total':14s} {m_tot:6.1f}           {M:6.1f}")

    rho = P * M / (R * T_GAS)
    mdot = m_tot * 1000 / 3600
    Q = mdot / rho
    A_in, A_th = math.pi * D_IN**2 / 4, math.pi * D_THROAT**2 / 4
    a = math.sqrt(gamma * R / M * T_GAS)
    print(f"\npseudo-gas: M = {M:.1f} kg/kmol, cp = {cp:.0f} J/(kg K), gamma = {gamma:.3f}")
    print(f"            mu = {mu:.2e} Pa s, k = {kappa:.3f} W/(m K), Pr = {mu * cp / kappa:.2f}")
    print(f"at {T_GAS - 273.15:.0f} C, {P / 1e3:.1f} kPa: rho = {rho:.3f} kg/m3, mdot = {mdot:.2f} kg/s, "
          f"Q = {Q:.2f} m3/s ({Q * 3600:.0f} m3/h)")
    print(f"velocity: inlet {Q / A_in:.1f} m/s, throat {Q / A_th:.1f} m/s, Mach throat {Q / A_th / a:.2f}")
    print(f"Re throat = {rho * Q / A_th * D_THROAT / mu:.2e}")

    m_oil = OIL_TPH * 1000 / 3600
    q_oil = m_oil / OIL["rho"]
    v_jet = q_oil / (math.pi * D_ORIFICE**2 / 4)
    print(f"\noil: {m_oil:.2f} kg/s, {q_oil * 1000:.1f} l/s; L/G = {q_oil / Q * 1000:.2f} l/m3, "
          f"mass ratio {m_oil / mdot:.2f}")
    print(f"orifice d = {D_ORIFICE * 1000:.0f} mm: jet {v_jet:.1f} m/s, dp_nozzle ~ {0.5 * OIL['rho'] * v_jet**2 / 0.8**2 / 1e5:.1f} bar")
    T_mix = (mdot * cp * T_GAS + m_oil * OIL["cp"] * T_OIL) / (mdot * cp + m_oil * OIL["cp"])
    print(f"adiabatic mixing, no phase change: T = {T_mix - 273.15:.0f} C (upper bound of cooling)")
    v_th = Q / A_th
    dp_gas = 0.5 * rho * v_th**2
    dp_liq = OIL["rho"] * v_th**2 * q_oil / Q   # Calvert-type acceleration loss
    print(f"pressure drop estimate: ~{(dp_liq + 0.3 * dp_gas) / 1e3:.1f} kPa "
          f"(dyn. head throat {dp_gas / 1e3:.1f} kPa, liquid acceleration {dp_liq / 1e3:.1f} kPa)")
    # Nukiyama-Tanasawa, CGS units: sigma dyn/cm, rho g/cm3, mu poise, L/G in l/m3 (=1000*QL/QG)
    sig, rl, mul = OIL["sigma"] * 1e3, OIL["rho"] * 1e-3, OIL["mu"] * 10
    d32 = 585 / v_th * math.sqrt(sig / rl) + 597 * (mul / math.sqrt(sig * rl))**0.45 * (1000 * q_oil / Q)**1.5
    print(f"Nukiyama-Tanasawa d32 ~ {d32:.0f} um (rough estimate for the droplet size distribution)")

if __name__ == "__main__":
    main()
