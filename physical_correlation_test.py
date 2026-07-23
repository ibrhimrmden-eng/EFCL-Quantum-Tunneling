import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.stats import pearsonr, spearmanr
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T
from lambda_star_extrapolation_test import find_lambda_star

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0

def V_std_scalar(x):
    return V_barrier_smooth(np.array([x]), V0=V0, a=a, w=w)[0]

def turning_point(E):
    """x2>0 where V_std(x2) = E (V_std is symmetric, monotonically decreasing for x>0)."""
    if E >= V_std_scalar(0.0):
        return np.nan
    return brentq(lambda x: V_std_scalar(x) - E, 0.0, 6.0, xtol=1e-9)

def S0_L_Vbar(E):
    x2 = turning_point(E)
    if np.isnan(x2):
        return np.nan, np.nan, np.nan
    x1 = -x2
    L = x2 - x1
    S0, _ = quad(lambda x: np.sqrt(max(0.0, 2 * (V_std_scalar(x) - E))), x1, x2, limit=200)
    Vint, _ = quad(V_std_scalar, x1, x2, limit=200)
    Vbar = Vint / L
    return S0, L, Vbar

# ------------------------------------------------------------
# Energy sample: only E < V0 (turning points must exist for the bare barrier)
# ------------------------------------------------------------
E_list = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99]

print(f"{'E':>6} {'lambda*':>10} {'S0(E)':>10} {'L(E)':>10} {'Vbar(E)':>10}")
data = []
for E in E_list:
    ls, T0 = find_lambda_star(E)
    S0, L, Vbar = S0_L_Vbar(E)
    print(f"{E:6.2f} {ls:10.5f} {S0:10.5f} {L:10.5f} {Vbar:10.5f}")
    data.append((E, ls, S0, L, Vbar))

data = np.array(data)
E_arr, lam_star, S0_arr, L_arr, Vbar_arr = data.T

print()
print("=" * 90)
print("Correlation of lambda*(E) with each candidate physical quantity")
print("=" * 90)

def report_corr(name, x, y):
    pear, p_p = pearsonr(x, y)
    spear, p_s = spearmanr(x, y)
    b, aa = np.polyfit(x, y, 1)
    pred = aa + b * x
    rmse = np.sqrt(np.mean((pred - y) ** 2))
    print(f"{name:>10}: Pearson r={pear:+.4f} (p={p_p:.2e})   Spearman rho={spear:+.4f} (p={p_s:.2e})   "
          f"linear-fit RMSE={rmse:.5f}")
    return pear, spear, rmse

print("\n-- direct correlations with lambda*(E) --")
res_S0 = report_corr("S0", S0_arr, lam_star)
res_L = report_corr("L", L_arr, lam_star)
res_Vbar = report_corr("Vbar", Vbar_arr, lam_star)

print("\n-- monotonic transforms tested --")
report_corr("1/S0", 1.0 / S0_arr, lam_star)
report_corr("sqrt(S0)", np.sqrt(S0_arr), lam_star)
report_corr("1/L", 1.0 / L_arr, lam_star)
report_corr("Vbar^2", Vbar_arr ** 2, lam_star)
res_E = report_corr("E (for reference)", E_arr, lam_star)

print()
print("=" * 90)
print("Reference: how good was the quadratic-in-E fit on this same E-range for comparison?")
print("=" * 90)
c0, c1, c2 = 0.2178, 0.0600, 0.0766
pred_quadE = c0 + c1 * E_arr + c2 * E_arr ** 2
rmse_quadE = np.sqrt(np.mean((pred_quadE - lam_star) ** 2))
print(f"quadratic-in-E RMSE on this E-range: {rmse_quadE:.5f}")

# --- CSV export (added for figure/table reproducibility; uses the actual
#     computed return values above, not transcribed numbers) ---
import csv
with open('data_csv/table4_semiclassical.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['quantity', 'pearson_r', 'spearman_rho', 'rmse'])
    writer.writerow(['E_raw', res_E[0], res_E[1], res_E[2]])
    writer.writerow(['S0', res_S0[0], res_S0[1], res_S0[2]])
    writer.writerow(['L', res_L[0], res_L[1], res_L[2]])
    writer.writerow(['Vbar', res_Vbar[0], res_Vbar[1], res_Vbar[2]])
print("Saved: data_csv/table4_semiclassical.csv")
