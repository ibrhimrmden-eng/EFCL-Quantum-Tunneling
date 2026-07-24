import numpy as np
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

def V_efcl_lam(x, lam1):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1 * env_field(x, t=0.0, sigma=sigma)

E_vals = [0.2, 0.5, 0.8, 0.95, 1.0, 1.1]
lam_grid = np.array([0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.65,0.80,1.00,1.3,1.6,2.0])
thresholds = [-0.5, -1.0, -1.5, -2.0]

results = {}
for E in E_vals:
    T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
    Rvals = []
    for lam in lam_grid:
        Te, _ = logderiv_T(lambda x, **kw: V_efcl_lam(x, lam1=lam), E, x_L=xL, x_R=xR)
        Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
    Rvals = np.array(Rvals)
    lam_star_row = []
    for t in thresholds:
        if np.nanmin(Rvals) < t < np.nanmax(Rvals):
            ls = np.interp(t, Rvals[::-1], lam_grid[::-1])
        else:
            ls = np.nan
        lam_star_row.append(ls)
    results[E] = lam_star_row
    print(f"E={E:.2f}  T0={T0:.4e}  lambda* at R={thresholds}: {['%.5f'%v if not np.isnan(v) else 'N/A' for v in lam_star_row]}")

print()
print("Ratios across thresholds (should be E-independent if single-variable collapse holds):")
for i,t in enumerate(thresholds):
    if t==-1.0: continue
    ratios=[]
    for E in E_vals:
        l_t = results[E][i]
        l_ref = results[E][thresholds.index(-1.0)]
        if not (np.isnan(l_t) or np.isnan(l_ref)):
            ratios.append(l_t/l_ref)
    ratios=np.array(ratios)
    if len(ratios)>1:
        print(f"R={t} vs R=-1.0: ratios={ratios}  mean={ratios.mean():.5f} CV={ratios.std()/ratios.mean():.4f}")

# --- CSV export (added for figure/table reproducibility; does not alter computation above) ---
import csv
with open('data_csv/table5_threshold_invariance.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['E'] + [f'lambda_star_R={t}' for t in thresholds])
    for E in E_vals:
        writer.writerow([E] + results[E])
print("Saved: data_csv/table5_threshold_invariance.csv")
