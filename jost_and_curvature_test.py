import numpy as np
from scipy.stats import pearsonr, spearmanr
from logderiv_ref import V_barrier_smooth, logderiv_T

V0, a, w = 1.0, 1.0, 0.12
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

E_list = np.array([0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99])
lam_star = np.array([0.20330,0.20617,0.20903,0.21609,0.22991,0.24332,0.25640,0.26940,0.28266,0.29675,0.31243,0.33069,0.34118,0.34799,0.35036])
delta = np.array([1.52801,1.49733,1.47680,1.44090,1.39623,1.36912,1.35310,1.34581,1.34620,1.35371,1.36799,1.38869,1.40131,1.40956,1.41241])

T0_vals = []
for E in E_list:
    T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
    T0_vals.append(T0)
T0_vals = np.array(T0_vals)
F_mag = 1.0/np.sqrt(T0_vals)   # |F(k)| in the convention T=1/|F|^2

def report(name, x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    pear, pp = pearsonr(x, y)
    spear, ps = spearmanr(x, y)
    b, aa = np.polyfit(x, y, 1)
    pred = aa + b*x
    rmse = np.sqrt(np.mean((pred-y)**2))
    print(f"{name:>22}: Pearson r={pear:+.4f} (p={pp:.2e})  Spearman rho={spear:+.4f}  RMSE={rmse:.5f}")
    return pear, spear, rmse

print("=== |F(k)| ~ 1/sqrt(T0(E))  (simplest Jost-magnitude proxy) vs lambda*(E) ===")
for E,f,l in zip(E_list, F_mag, lam_star):
    print(f"  E={E:.2f}  T0={T0_vals[list(E_list).index(E)]:.4e}  |F|={f:.4f}  lambda*={l:.5f}")
res_Fmag = report("|F(k)|", F_mag, lam_star)
res_lnFmag = report("ln|F(k)| (~ -0.5 lnT0)", np.log(F_mag), lam_star)
report("T0(E) itself", T0_vals, lam_star)

# --- CSV export (added for figure/table reproducibility; uses the actual
#     computed return values above, not transcribed numbers) ---
import csv
with open('data_csv/table4_jost.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['quantity', 'pearson_r', 'spearman_rho', 'rmse'])
    writer.writerow(['F_mag', res_Fmag[0], res_Fmag[1], res_Fmag[2]])
    writer.writerow(['ln_F_mag', res_lnFmag[0], res_lnFmag[1], res_lnFmag[2]])
print("Saved: data_csv/table4_jost.csv")

print()
print("=== curvature test: d2(delta)/dE2 vs dlambda*/dE  (feature-location comparison) ===")
d_delta_dE = np.gradient(delta, E_list)
d2_delta_dE2 = np.gradient(d_delta_dE, E_list)
d_lamstar_dE = np.gradient(lam_star, E_list)

print(f"{'E':>6} {'d2delta/dE2':>14} {'dlambda*/dE':>14}")
for E, d2, dl in zip(E_list, d2_delta_dE2, d_lamstar_dE):
    print(f"{E:6.2f} {d2:14.5f} {dl:14.5f}")

report("d2(delta)/dE2", d2_delta_dE2, d_lamstar_dE)
