import numpy as np
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

V0, a, sigma = 1.0, 1.0, 1.2
xL, xR = -8.0, 8.0

def make_Vfuncs(w):
    def V_std(x, **kw):
        return V_barrier_smooth(x, V0=V0, a=a, w=w)
    def V_efcl_lam(x, lam1):
        return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1 * env_field(x, t=0.0, sigma=sigma)
    return V_std, V_efcl_lam

E_list = [0.05, 0.20, 0.40, 0.60, 0.80, 0.95, 0.99]
lam_grid = np.array([0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.65,0.80,1.00])

results = {}
for w in [0.08, 0.12, 0.20]:
    V_std, V_efcl_lam = make_Vfuncs(w)
    lam_star_row = []
    for E in E_list:
        T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
        Rvals = []
        for lam in lam_grid:
            Te, _ = logderiv_T(lambda x, **kw: V_efcl_lam(x, lam1=lam), E, x_L=xL, x_R=xR)
            Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
        Rvals = np.array(Rvals)
        if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
            ls = np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
        else:
            ls = np.nan
        lam_star_row.append(ls)
    results[w] = np.array(lam_star_row)
    print(f"w={w}: lambda*(E) = {['%.5f'%v if not np.isnan(v) else 'NA' for v in lam_star_row]}")

print()
print("=" * 90)
print(f"{'E':>6} {'lam*(w=0.08)':>14} {'lam*(w=0.12)':>14} {'lam*(w=0.20)':>14} {'ratio(0.08/0.12)':>18} {'ratio(0.20/0.12)':>18}")
ratios_08 = []
ratios_20 = []
for i, E in enumerate(E_list):
    l08, l12, l20 = results[0.08][i], results[0.12][i], results[0.20][i]
    r1 = l08/l12
    r2 = l20/l12
    ratios_08.append(r1); ratios_20.append(r2)
    print(f"{E:6.2f} {l08:14.5f} {l12:14.5f} {l20:14.5f} {r1:18.5f} {r2:18.5f}")

ratios_08 = np.array(ratios_08); ratios_20 = np.array(ratios_20)
print()
print(f"ratio(w=0.08/w=0.12): mean={ratios_08.mean():.4f} std={ratios_08.std():.4f} CV={ratios_08.std()/ratios_08.mean():.4f}")
print(f"ratio(w=0.20/w=0.12): mean={ratios_20.mean():.4f} std={ratios_20.std():.4f} CV={ratios_20.std()/ratios_20.mean():.4f}")

# --- CSV export (added for figure/table reproducibility; does not alter computation above) ---
import csv
with open('data_csv/table6_barrier_width.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['E'] + [f'lambda_star_w={wv}' for wv in [0.08,0.12,0.20]])
    for i, E in enumerate(E_list):
        writer.writerow([E, results[0.08][i], results[0.12][i], results[0.20][i]])
print("Saved: data_csv/table6_barrier_width.csv")
