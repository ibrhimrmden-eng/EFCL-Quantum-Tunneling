import numpy as np
from logderiv_ref import V_barrier_smooth, logderiv_T

V0, a, w = 1.0, 1.0, 0.12
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

def env_field_sigma(x, sigma):
    return np.exp(-(x/sigma)**2)

def V_efcl_lam(x, lam1, sigma):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1*env_field_sigma(x, sigma)

lam_grid = np.array([0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.65,0.80,1.00])
E_list = [0.2, 0.5, 0.8, 0.95]

results = {}
for sigma in [1.0, 1.2, 1.4]:
    lam_star_row = []
    for E in E_list:
        T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
        Rvals = []
        for lam in lam_grid:
            Te, _ = logderiv_T(lambda x, **kw: V_efcl_lam(x, lam1=lam, sigma=sigma), E, x_L=xL, x_R=xR)
            Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
        Rvals = np.array(Rvals)
        if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
            ls = np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
        else:
            ls = np.nan
        lam_star_row.append(ls)
    results[sigma] = np.array(lam_star_row)
    print(f"sigma={sigma}: lambda*(E) = {['%.5f'%v if not np.isnan(v) else 'NA' for v in lam_star_row]}")

print()
print(f"{'E':>6} {'lam*(s=1.0)':>13} {'lam*(s=1.2)':>13} {'lam*(s=1.4)':>13} {'ratio(1.0/1.2)':>16} {'ratio(1.4/1.2)':>16}")
r10=[]; r14=[]
for i, E in enumerate(E_list):
    l10, l12, l14 = results[1.0][i], results[1.2][i], results[1.4][i]
    a1 = l10/l12; a2 = l14/l12
    r10.append(a1); r14.append(a2)
    print(f"{E:6.2f} {l10:13.5f} {l12:13.5f} {l14:13.5f} {a1:16.5f} {a2:16.5f}")
r10=np.array(r10); r14=np.array(r14)
print()
print(f"ratio(sigma=1.0/1.2): mean={r10.mean():.4f} CV={r10.std()/r10.mean():.4f}")
print(f"ratio(sigma=1.4/1.2): mean={r14.mean():.4f} CV={r14.std()/r14.mean():.4f}")

# Also test envelope CENTER shift (small offset x0) as a second physical-noise check
print()
print("=== envelope center shift test (x0 = -0.3, 0, +0.3), sigma=1.2 fixed ===")
def V_efcl_shift(x, lam1, x0):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1*np.exp(-((x-x0)/1.2)**2)

results_shift = {}
for x0 in [-0.3, 0.0, 0.3]:
    lam_star_row = []
    for E in E_list:
        T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
        Rvals = []
        for lam in lam_grid:
            Te, _ = logderiv_T(lambda x, **kw: V_efcl_shift(x, lam1=lam, x0=x0), E, x_L=xL, x_R=xR)
            Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
        Rvals = np.array(Rvals)
        if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
            ls = np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
        else:
            ls = np.nan
        lam_star_row.append(ls)
    results_shift[x0] = np.array(lam_star_row)
    print(f"x0={x0}: lambda*(E) = {['%.5f'%v if not np.isnan(v) else 'NA' for v in lam_star_row]}")

r_m03=[]; r_p03=[]
for i, E in enumerate(E_list):
    lm, l0, lp = results_shift[-0.3][i], results_shift[0.0][i], results_shift[0.3][i]
    r_m03.append(lm/l0); r_p03.append(lp/l0)
r_m03=np.array(r_m03); r_p03=np.array(r_p03)
print(f"ratio(x0=-0.3/x0=0): mean={r_m03.mean():.4f} CV={r_m03.std()/r_m03.mean():.4f}")
print(f"ratio(x0=+0.3/x0=0): mean={r_p03.mean():.4f} CV={r_p03.std()/r_p03.mean():.4f}")

# --- CSV export (added for figure/table reproducibility; does not alter computation above) ---
import csv
with open('data_csv/table6_envelope_sigma.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['E'] + [f'lambda_star_sigma={s}' for s in [1.0,1.2,1.4]])
    for i, E in enumerate(E_list):
        writer.writerow([E, results[1.0][i], results[1.2][i], results[1.4][i]])
print("Saved: data_csv/table6_envelope_sigma.csv")
