import numpy as np
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

a, sigma, w = 1.0, 1.2, 0.12
xL, xR = -8.0, 8.0

def make_Vfuncs(V0):
    def V_std(x, **kw):
        return V_barrier_smooth(x, V0=V0, a=a, w=w)
    def V_efcl_lam(x, lam1):
        return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1 * env_field(x, t=0.0, sigma=sigma)
    return V_std, V_efcl_lam

lam_grid = np.array([0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.65,0.80,1.00,1.2,1.5])

# Use FRACTION of V0 for E, so we compare the same relative tunneling regime across barrier heights
frac_list = [0.05, 0.20, 0.40, 0.60, 0.80, 0.95]

results = {}
for V0 in [0.8, 1.0, 1.2]:
    V_std, V_efcl_lam = make_Vfuncs(V0)
    lam_star_row = []
    for frac in frac_list:
        E = frac * V0
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
    results[V0] = np.array(lam_star_row)
    print(f"V0={V0}: lambda*(E=frac*V0) = {['%.5f'%v if not np.isnan(v) else 'NA' for v in lam_star_row]}")

print()
print(f"{'frac':>6} {'lam*(V0=0.8)':>14} {'lam*(V0=1.0)':>14} {'lam*(V0=1.2)':>14} {'ratio(0.8/1.0)':>16} {'ratio(1.2/1.0)':>16}")
r08=[]; r12=[]
for i, frac in enumerate(frac_list):
    l08, l10, l12v = results[0.8][i], results[1.0][i], results[1.2][i]
    a1 = l08/l10; a2 = l12v/l10
    r08.append(a1); r12.append(a2)
    print(f"{frac:6.2f} {l08:14.5f} {l10:14.5f} {l12v:14.5f} {a1:16.5f} {a2:16.5f}")
r08=np.array(r08); r12=np.array(r12)
print()
print(f"ratio(V0=0.8/V0=1.0): mean={r08.mean():.4f} CV={r08.std()/r08.mean():.4f}")
print(f"ratio(V0=1.2/V0=1.0): mean={r12.mean():.4f} CV={r12.std()/r12.mean():.4f}")

# --- CSV export (added for figure/table reproducibility; does not alter computation above) ---
import csv
with open('data_csv/table6_barrier_V0.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['frac_E_over_V0'] + [f'lambda_star_V0={v}' for v in [0.8,1.0,1.2]])
    for i, frac in enumerate(frac_list):
        writer.writerow([frac, results[0.8][i], results[1.0][i], results[1.2][i]])
print("Saved: data_csv/table6_barrier_V0.csv")
