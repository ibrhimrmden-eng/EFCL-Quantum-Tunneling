import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logderiv_ref import V_barrier_smooth, env_field, logderiv_T

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)
def efcl_potential_localized(x, lam1=0.35, lam2=0.10, sigma=1.2):
    envelope = np.exp(-(x/sigma)**2)
    env = env_field(x, t=0.0, sigma=sigma)
    nu2 = (1.0+0.2*np.cos(0.0))**2
    return lam1*env + lam2*nu2*envelope
def V_efcl(x, **kw):
    return V_std(x) + efcl_potential_localized(x)

def wkb_transmission(x, Veff, E):
    forbidden = Veff > E
    if not np.any(forbidden):
        return 1.0
    kappa = np.sqrt(np.maximum(0.0, 2.0*(Veff-E)))
    integral = np.trapezoid(kappa[forbidden], x[forbidden])
    return float(np.clip(np.exp(-2.0*integral), 0.0, 1.0))

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

# ---------- Figure 1: Barrier shape before/after EFCL ----------
x = np.linspace(-5, 5, 2000)
Vs = V_std(x)
Ve = V_efcl(x)
fig, ax = plt.subplots(figsize=(6,4.2))
ax.plot(x, Vs, label=r'$V_{std}(x)$ (bare barrier)', lw=2, color='#1f77b4')
ax.plot(x, Ve, label=r'$V_{eff}(x)=V_{std}+\Phi_{EFCL}$', lw=2, color='#d62728', ls='--')
ax.set_xlabel('x')
ax.set_ylabel('Potential V(x)')
ax.set_title('Figure 1: Barrier profile before and after localized EFCL correction')
ax.legend(frameon=False)
ax.set_xlim(-5,5)
plt.tight_layout()
plt.savefig('figures/Figure1_barrier_profile.png', dpi=600)
plt.close()
print("Figure 1 done")

# ---------- Figure 2: Exact (Riccati) vs WKB, T(E) ----------
E_arr = np.linspace(0.05, 1.7, 60)
T_exact = []
T_wkb = []
xg = np.linspace(xL, xR, 6000)
Vsg = V_std(xg)
for E in E_arr:
    T0, _ = logderiv_T(V_std, E, x_L=xL, x_R=xR)
    T_exact.append(T0)
    T_wkb.append(wkb_transmission(xg, Vsg, E))
T_exact = np.array(T_exact); T_wkb = np.array(T_wkb)

fig, ax = plt.subplots(figsize=(6,4.2))
ax.plot(E_arr, T_exact, label='Exact (log-derivative)', color='#1f77b4', lw=2)
ax.plot(E_arr, T_wkb, label='WKB (primitive exponent)', color='#ff7f0e', lw=2, ls='--')
ax.axvline(V0, color='gray', lw=1, ls=':', label=r'$V_0$ (barrier top)')
ax.set_yscale('log')
ax.set_xlabel('Energy E')
ax.set_ylabel('Transmission T(E)')
ax.set_title('Figure 2: Exact transmission vs. primitive WKB')
ax.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig('figures/Figure2_exact_vs_WKB.png', dpi=600)
plt.close()
print("Figure 2 done")

# ---------- Figure 3: Grid convergence (hardcoded from validated session results) ----------
N_vals = np.array([500, 1000, 2000, 4000, 8000, 16000])
ratio_vals = np.array([0.09295736, 0.09331310, 0.09319245, 0.09311626, 0.09305363, 0.09306112])
maxPhi_vals = np.array([0.493976, 0.493994, 0.493998, 0.494000, 0.494000, 0.494000])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4))
ax1.semilogx(N_vals, ratio_vals, 'o-', color='#2ca02c')
ax1.set_xlabel('Grid points N')
ax1.set_ylabel(r'$T_{EFCL}/T_{std}$ at E=0.9')
ax1.set_title('(a) Suppression ratio convergence')
ax1.ticklabel_format(axis='y', style='plain')

ax2.semilogx(N_vals, maxPhi_vals, 's-', color='#9467bd')
ax2.set_xlabel('Grid points N')
ax2.set_ylabel(r'max $\Phi_{EFCL}(x)$')
ax2.set_title('(b) Correction amplitude convergence')

fig.suptitle('Figure 3: Grid convergence after localization fix')
plt.tight_layout()
plt.savefig('figures/Figure3_grid_convergence.png', dpi=600)
plt.close()
print("Figure 3 done")
