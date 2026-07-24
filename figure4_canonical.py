"""
CANONICAL SCRIPT for Figure 4 (scaling collapse).
This is the single source of truth for Figure 4's data. Any change to
Figure 4 must be made by editing and re-running this script, then
regenerating the figure from its CSV output -- never by editing the
plot directly.
"""
import numpy as np
import csv
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)
def V_efcl_lam(x, lam1):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1*env_field(x, t=0.0, sigma=sigma)

E_list = [0.2, 0.5, 0.8, 0.95, 1.05]
lam_grid = np.array([0.02,0.05,0.08,0.12,0.16,0.20,0.30,0.40,0.55,0.70])

rows = []
for E in E_list:
    T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
    Rvals = []
    for lam in lam_grid:
        Te, _ = logderiv_T(lambda x, **kw: V_efcl_lam(x, lam1=lam), E, x_L=xL, x_R=xR)
        Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
    Rvals = np.array(Rvals)

    if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
        lam_star = np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
    else:
        lam_star = np.nan

    for lam, R in zip(lam_grid, Rvals):
        rows.append([E, lam, R, lam_star, lam/lam_star if lam_star and not np.isnan(lam_star) else np.nan])
    print(f"E={E}: lambda*={lam_star:.5f}")

with open('data_csv/figure4_scaling_collapse.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['E', 'lambda', 'R', 'lambda_star', 'lambda_over_lambda_star'])
    writer.writerows(rows)
print("Saved: data_csv/figure4_scaling_collapse.csv")
