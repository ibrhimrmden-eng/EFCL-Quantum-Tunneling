import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr, spearmanr
from logderiv_ref import V_barrier_smooth, env_field

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0
hbar, m = 1.0, 1.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)

def logderiv_full(V_func, E, x_L=-8.0, x_R=8.0, **Vkwargs):
    """Same as logderiv_T but also returns the complex reflection amplitude r."""
    V_R = V_func(np.array([x_R]), **Vkwargs)[0]
    V_L = V_func(np.array([x_L]), **Vkwargs)[0]
    k_out = np.sqrt(2*m*complex(E-V_R))/hbar
    k_in  = np.sqrt(2*m*complex(E-V_L))/hbar
    if k_out == 0: k_out = 1e-9
    if k_in == 0: k_in = 1e-9

    def rhs(x, yc):
        y = yc[0] + 1j*yc[1]
        Vx = V_func(np.array([x]), **Vkwargs)[0]
        dy = (2*m/hbar**2)*(Vx-E) - y**2
        return [dy.real, dy.imag]

    y0 = [np.real(1j*k_out), np.imag(1j*k_out)]
    sol = solve_ivp(rhs, [x_R, x_L], y0, method="RK45", rtol=1e-10, atol=1e-12, max_step=(x_R-x_L)/4000)
    y_L = sol.y[0,-1] + 1j*sol.y[1,-1]
    r = np.exp(2j*k_in*x_L) * (1j*k_in - y_L) / (1j*k_in + y_L)
    R = np.abs(r)**2
    T = 1.0 - R
    return T, R, r

# ------------------------------------------------------------
# Background (bare barrier) scattering phase delta_std(E) = 0.5*arg(r_std)
# ------------------------------------------------------------
E_list = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99]
lam_star = np.array([0.20330,0.20617,0.20903,0.21609,0.22991,0.24332,0.25640,0.26940,0.28266,0.29675,0.31243,0.33069,0.34118,0.34799,0.35036])

deltas = []
for E in E_list:
    T, R, r = logderiv_full(V_std, E, x_L=xL, x_R=xR)
    delta = 0.5*np.angle(r)
    deltas.append(delta)
deltas = np.array(deltas)

# unwrap phase for smooth dDelta/dE
deltas_unwrapped = np.unwrap(deltas)
dDelta_dE = np.gradient(deltas_unwrapped, E_list)

print(f"{'E':>6} {'lambda*':>10} {'delta(E)':>10} {'delta_unwr':>12} {'dDelta/dE':>12}")
for E, ls, d, du, dd in zip(E_list, lam_star, deltas, deltas_unwrapped, dDelta_dE):
    print(f"{E:6.2f} {ls:10.5f} {d:10.5f} {du:12.5f} {dd:12.5f}")

def report(name, x, y):
    x = np.asarray(x); y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    pear, pp = pearsonr(x, y)
    spear, ps = spearmanr(x, y)
    b, aa = np.polyfit(x, y, 1)
    pred = aa + b*x
    rmse = np.sqrt(np.mean((pred-y)**2))
    print(f"{name:>18}: Pearson r={pear:+.4f} (p={pp:.2e})  Spearman rho={spear:+.4f}  RMSE(linear fit)={rmse:.5f}")
    return pear, spear, rmse

print()
print("=== correlations of lambda*(E) with scattering-phase-derived quantities ===")
res_delta = report("delta", deltas_unwrapped, lam_star)
res_tandelta = report("tan(delta)", np.tan(deltas), lam_star)
report("sin(delta)", np.sin(deltas), lam_star)
report("cos(delta)", np.cos(deltas), lam_star)
res_wigner = report("dDelta/dE (Wigner delay)", dDelta_dE, lam_star)

# --- CSV export (added for figure/table reproducibility; uses the actual
#     computed return values above, not transcribed numbers) ---
import csv
with open('data_csv/table4_phase_shift.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['quantity', 'pearson_r', 'spearman_rho', 'rmse'])
    writer.writerow(['delta', res_delta[0], res_delta[1], res_delta[2]])
    writer.writerow(['tan_delta', res_tandelta[0], res_tandelta[1], res_tandelta[2]])
    writer.writerow(['wigner_delay', res_wigner[0], res_wigner[1], res_wigner[2]])
print("Saved: data_csv/table4_phase_shift.csv")
