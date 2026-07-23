import numpy as np
from scipy.optimize import brentq
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

def V_efcl_lam(x, lam1):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1 * env_field(x, t=0.0, sigma=sigma)

def find_lambda_star(E, target_R=-1.0, lam_lo=1e-4, lam_hi=3.0):
    """Root-find the lambda where R(E,lambda) = target_R, expanding bracket if needed."""
    T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)

    def g(lam):
        Te, _ = logderiv_T(lambda x, **kw: V_efcl_lam(x, lam1=lam), E, x_L=xL, x_R=xR)
        if Te <= 0:
            return -50.0 - target_R  # very negative R, past the root
        return np.log(Te / T0) - target_R

    lo, hi = lam_lo, lam_hi
    glo, ghi = g(lo), g(hi)
    tries = 0
    while glo * ghi > 0 and tries < 10:
        hi *= 1.7
        ghi = g(hi)
        tries += 1
    if glo * ghi > 0:
        return np.nan, T0
    root = brentq(g, lo, hi, xtol=1e-6, rtol=1e-8)
    return root, T0

# ------------------------------------------------------------
# The quadratic model fitted PREVIOUSLY on E in [0.2, 1.1]  (do not refit)
# lambda*(E) = c0 + c1*E + c2*E^2
# ------------------------------------------------------------
c0, c1, c2 = 0.2178, 0.0600, 0.0766
def quad_model(E):
    return c0 + c1 * E + c2 * E**2

regions = {
    "Region 1 -- very low E (below training range, deep tunneling)":
        [0.01, 0.03, 0.05, 0.10],
    "Region 2 -- around barrier top (crossover region, mostly interpolation)":
        [0.90, 0.95, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05, 1.10],
    "Region 3 -- above barrier (genuine extrapolation, over-barrier regime)":
        [1.20, 1.40, 1.60, 1.80, 2.00],
}

all_results = []
for region_name, E_list in regions.items():
    print("=" * 100)
    print(region_name)
    print("=" * 100)
    print(f"{'E':>8} {'lambda*_exact':>16} {'T0(E)':>14} {'lambda*_quad_pred':>18} {'abs error':>12}")
    for E in E_list:
        try:
            ls_exact, T0 = find_lambda_star(E)
        except Exception as ex:
            ls_exact, T0 = np.nan, np.nan
        ls_pred = quad_model(E)
        err = abs(ls_pred - ls_exact) if not np.isnan(ls_exact) else np.nan
        print(f"{E:8.2f} {ls_exact:16.5f} {T0:14.4e} {ls_pred:18.5f} {err:12.5f}")
        all_results.append((region_name, E, ls_exact, ls_pred, err))

# ------------------------------------------------------------
# RMSE per region, and for TRUE extrapolation only (regions 1 and 3)
# ------------------------------------------------------------
print()
print("=" * 100)
print("RMSE summary")
print("=" * 100)
for region_name in regions:
    errs = [r[4] for r in all_results if r[0] == region_name and not np.isnan(r[4])]
    if errs:
        rmse = np.sqrt(np.mean(np.array(errs) ** 2))
        print(f"{region_name}: RMSE = {rmse:.5f}  (n={len(errs)})")

true_extrap_errs = [r[4] for r in all_results
                     if ("Region 1" in r[0] or "Region 3" in r[0]) and not np.isnan(r[4])]
if true_extrap_errs:
    rmse_extrap = np.sqrt(np.mean(np.array(true_extrap_errs) ** 2))
    print(f"\nTRUE OUT-OF-SAMPLE (Region 1 + Region 3 combined): RMSE = {rmse_extrap:.5f}  (n={len(true_extrap_errs)})")

# for reference, in-sample RMSE (training region E=0.2-1.1, from previous fit)
train_E = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
train_lam_star = np.array([0.22993, 0.24332, 0.25641, 0.26940, 0.28266, 0.29675, 0.31242, 0.33069, 0.35277, 0.38010])
train_pred = quad_model(train_E)
rmse_train = np.sqrt(np.mean((train_pred - train_lam_star) ** 2))
print(f"IN-SAMPLE (training region E=0.2-1.1): RMSE = {rmse_train:.5f}  (n={len(train_E)})  [reference baseline]")
