"""Droplet size of the oil spray from engineering correlations.

1. Sheet from the deflector disk: thickness h = Q / (pi D U) at the rim.
2. Secondary breakup by the gas: the largest stable droplet satisfies
   We = rho_g U_rel^2 d / sigma = We_cr,  We_cr = 12 (1 + 1.077 Oh^1.6)  (Brodkey),
   Oh = mu / sqrt(rho_l sigma d)  -> solved iteratively; d32 ~ 0.4...0.6 d_max.
3. Venturi throat atomisation, Nukiyama-Tanasawa (CGS units):
   d32[um] = 585/U * sqrt(sigma/rho) + 597 (mu/sqrt(sigma rho))^0.45 (1000 QL/QG)^1.5

Edit the inputs below and run:  python droplet_size.py
"""
import math

# ---- inputs
M_OIL = 32000 / 3600      # kg/s
RHO_L = 1150.0            # kg/m3 at 100 C (1.2 kg/l cold)
NU_L = 40e-6              # m2/s at 100 C (40 cSt)
SIGMA = 0.030             # N/m, assumed (typical heavy oil)
D_DISK, U_SHEET = 0.065, 25.0
RHO_G, Q_G, U_THROAT = 1.106, 19.09, 77.5


def d_max(u_rel):
    d = 12 * SIGMA / (RHO_G * u_rel**2)
    for _ in range(50):
        oh = RHO_L * NU_L / math.sqrt(RHO_L * SIGMA * d)
        d = 12 * (1 + 1.077 * oh**1.6) * SIGMA / (RHO_G * u_rel**2)
    return d, oh


def main():
    mu = RHO_L * NU_L
    q = M_OIL / RHO_L
    print(f"oil: {q * 1000:.1f} l/s, mu = {mu * 1000:.0f} mPa s, sigma = {SIGMA} N/m (assumed)")
    print(f"sheet at the disk rim: {q / (math.pi * D_DISK * U_SHEET) * 1000:.2f} mm thick")
    for u in (50, 60, 70):
        d, oh = d_max(u)
        print(f"secondary breakup, U_rel = {u} m/s: d_max = {d * 1e6:4.0f} um (Oh = {oh:.2f}) "
              f"-> d32 ~ {0.4 * d * 1e6:.0f}...{0.6 * d * 1e6:.0f} um")
    s, r, m = SIGMA * 1e3, RHO_L * 1e-3, mu * 10
    d32 = 585 / U_THROAT * math.sqrt(s / r) + 597 * (m / math.sqrt(s * r))**0.45 * (1000 * q / Q_G)**1.5
    print(f"Nukiyama-Tanasawa (throat {U_THROAT} m/s, L/G {1000 * q / Q_G:.2f} l/m3): d32 = {d32:.0f} um")
    # Rosin-Rammler d (63 % size) giving this d32 for n = 2.5: d32 = d / Gamma(1 - 1/n)
    print(f"Rosin-Rammler for cloudProperties (n = 2.5): d = {d32 * math.gamma(1 - 1 / 2.5):.0f} um")


if __name__ == "__main__":
    main()
