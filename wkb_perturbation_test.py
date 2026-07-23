import numpy as np
from logderiv_ref import (V_barrier_smooth, env_field, frequency_term,
                           efcl_potential_localized, logderiv_T)

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0
hbar, m = 1.0, 1.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

def V_efcl(x, t=0.0, lam1=0.0, lam2=0.0):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + efcl_potential_localized(x, t=t, lam1=lam1, lam2=lam2, sigma=sigma)


def wkb_perturbation_coeffs(E, Phi_func, x_grid):
    Vstd = V_std(x_grid)
    forbidden = Vstd > E
    A = 2 * m * (Vstd[forbidden] - E) / hbar**2
    Phi = Phi_func(x_grid[forbidden])
    B = 2 * m * Phi / hbar**2
    A_min = np.min(A) if len(A) else np.nan
    I1 = np.trapezoid(B / np.sqrt(A), x_grid[forbidden])
    I2 = np.trapezoid(B**2 / (4 * A**1.5), x_grid[forbidden])
    return I1, I2, A_min


if __name__ == "__main__":
    x_grid = np.linspace(xL, xR, 20000)

    for E_test in [0.9, 0.3]:
        print("=" * 78)
        print(f"E = {E_test}  (lam1 sweep, lam2=0, Phi(x) = env_field(x,0))")
        print("=" * 78)

        Ts_ref, _ = logderiv_T(V_std, E_test, x_L=xL, x_R=xR)

        def Phi_env(x):
            return env_field(x, t=0.0, sigma=sigma)

        I1, I2, A_min = wkb_perturbation_coeffs(E_test, Phi_env, x_grid)
        print(f"predicted (WKB pert.): I1={I1:.5f}  I2={I2:.5f}   min(A) in forbidden region={A_min:.6f}")
        print(f"predicted ln(ratio) ~= -I1*lam + I2*lam^2")
        print()

        lam_vals = np.array([0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20])
        ln_ratios_exact = []
        for lam1 in lam_vals:
            Te, _ = logderiv_T(lambda x, **kw: V_efcl(x, t=0.0, lam1=lam1, lam2=0.0), E_test, x_L=xL, x_R=xR)
            ln_ratios_exact.append(np.log(Te / Ts_ref))
        ln_ratios_exact = np.array(ln_ratios_exact)

        ln_ratios_pred = -I1 * lam_vals + I2 * lam_vals**2

        print(f"{'lambda1':>8} {'ln(ratio)_exact':>18} {'ln(ratio)_WKBpred':>20} {'abs diff':>12}")
        for i, lam1 in enumerate(lam_vals):
            print(f"{lam1:8.3f} {ln_ratios_exact[i]:18.6f} {ln_ratios_pred[i]:20.6f} {abs(ln_ratios_exact[i]-ln_ratios_pred[i]):12.6f}")

        c2, c1, c0 = np.polyfit(lam_vals, ln_ratios_exact, 2)
        print()
        print(f"direct quadratic fit to EXACT data: ln(ratio) ~= {c1:.5f}*lam + {c2:.5f}*lam^2  (c0={c0:.2e})")
        print(f"WKB-predicted coefficients:                      {-I1:.5f}*lam + {I2:.5f}*lam^2")
        print()
