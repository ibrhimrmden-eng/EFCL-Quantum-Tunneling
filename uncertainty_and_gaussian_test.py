import numpy as np
from scipy.integrate import solve_ivp
from logderiv_ref import V_barrier_smooth, env_field

V0, a, w, sigma = 1.0, 1.0, 0.12, 1.2
xL, xR = -8.0, 8.0
hbar, m = 1.0, 1.0

def V_std(x, **kw):
    return V_barrier_smooth(x, V0=V0, a=a, w=w)
def V_efcl_lam(x, lam1):
    return V_barrier_smooth(x, V0=V0, a=a, w=w) + lam1 * env_field(x, t=0.0, sigma=sigma)

def logderiv_T_custom(V_func, E, x_L=-8.0, x_R=8.0, rtol=1e-10, atol=1e-12, max_step_div=4000, **Vkwargs):
    V_R = V_func(np.array([x_R]), **Vkwargs)[0]
    V_L = V_func(np.array([x_L]), **Vkwargs)[0]
    k_out = np.sqrt(2*m*complex(E-V_R))/hbar
    k_in  = np.sqrt(2*m*complex(E-V_L))/hbar
    if k_out == 0: k_out = 1e-9
    if k_in == 0: k_in = 1e-9
    def rhs(x, yc):
        y = yc[0]+1j*yc[1]
        Vx = V_func(np.array([x]), **Vkwargs)[0]
        dy = (2*m/hbar**2)*(Vx-E) - y**2
        return [dy.real, dy.imag]
    y0=[np.real(1j*k_out), np.imag(1j*k_out)]
    sol = solve_ivp(rhs, [x_R,x_L], y0, method="RK45", rtol=rtol, atol=atol, max_step=(x_R-x_L)/max_step_div)
    y_L = sol.y[0,-1]+1j*sol.y[1,-1]
    r = np.exp(2j*k_in*x_L)*(1j*k_in-y_L)/(1j*k_in+y_L)
    T = 1.0 - np.abs(r)**2
    return float(T)

def find_lambda_star_custom(E, lam_grid, **solver_kw):
    T0 = logderiv_T_custom(V_std, E, x_L=xL, x_R=xR, **solver_kw)
    Rvals=[]
    for lam in lam_grid:
        Te = logderiv_T_custom(lambda x,**kw: V_efcl_lam(x,lam1=lam), E, x_L=xL, x_R=xR, **solver_kw)
        Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
    Rvals=np.array(Rvals)
    if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
        return np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
    return np.nan

lam_grid = np.array([0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50])
E_test = [0.3, 0.7, 0.95]

print("=== NUMERICAL UNCERTAINTY / SENSITIVITY TEST ===")
settings = {
    "default (rtol=1e-10,atol=1e-12,steps=4000)": dict(rtol=1e-10, atol=1e-12, max_step_div=4000),
    "looser (rtol=1e-8,atol=1e-10,steps=1000)":   dict(rtol=1e-8,  atol=1e-10, max_step_div=1000),
    "tighter (rtol=1e-12,atol=1e-14,steps=8000)": dict(rtol=1e-12, atol=1e-14, max_step_div=8000),
}
results = {}
for name, kw in settings.items():
    row = [find_lambda_star_custom(E, lam_grid, **kw) for E in E_test]
    results[name] = row
    print(f"{name}: lambda* = {['%.6f'%v for v in row]}")

ref = np.array(results["default (rtol=1e-10,atol=1e-12,steps=4000)"])
for name in results:
    if name.startswith("default"): continue
    diffpct = 100*np.abs(np.array(results[name])-ref)/ref
    print(f"  vs default, % diff for '{name}': {diffpct}")

print()
print("=== DIFFERENT BARRIER FAMILY: Gaussian barrier V(x) = V0*exp(-x^2/(2 s^2)) ===")
s_gauss = 0.6  # chosen so width is roughly comparable to tanh barrier
def V_gauss(x, **kw):
    return V0*np.exp(-(x**2)/(2*s_gauss**2))
def V_gauss_efcl(x, lam1):
    return V_gauss(x) + lam1*env_field(x, t=0.0, sigma=sigma)

def logderiv_T_generic(V_func, E, x_L=-8.0, x_R=8.0, **Vkwargs):
    return logderiv_T_custom(V_func, E, x_L=x_L, x_R=x_R, rtol=1e-10, atol=1e-12, max_step_div=4000, **Vkwargs)

def find_lstar_generic(V_std_f, V_efcl_f, E, lam_grid):
    T0 = logderiv_T_generic(V_std_f, E, x_L=xL, x_R=xR)
    Rvals=[]
    for lam in lam_grid:
        Te = logderiv_T_generic(lambda x,**kw: V_efcl_f(x,lam1=lam), E, x_L=xL, x_R=xR)
        Rvals.append(np.log(Te/T0) if Te>0 else np.nan)
    Rvals=np.array(Rvals)
    if np.nanmin(Rvals) < -1 < np.nanmax(Rvals):
        return np.interp(-1.0, Rvals[::-1], lam_grid[::-1])
    return np.nan

E_gauss = [0.2, 0.5, 0.8]
lam_star_gauss = [find_lstar_generic(V_gauss, V_gauss_efcl, E, lam_grid) for E in E_gauss]
print(f"Gaussian barrier lambda*(E): E={E_gauss} -> lambda*={lam_star_gauss}")
# compare monotonic trend qualitative match with tanh-barrier lambda*
lam_star_tanh_ref = [find_lambda_star_custom(E, lam_grid, rtol=1e-10, atol=1e-12, max_step_div=4000) for E in E_gauss]
print(f"tanh barrier (same E, V0=1) lambda*(E) for comparison: {lam_star_tanh_ref}")
