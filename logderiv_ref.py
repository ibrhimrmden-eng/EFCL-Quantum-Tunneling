import numpy as np
from scipy.integrate import solve_ivp

hbar = 1.0
m = 1.0

def V_barrier_smooth(x, V0=1.0, a=1.0, w=0.12):
    return 0.5 * V0 * (np.tanh((x + a) / w) - np.tanh((x - a) / w))

def env_field(x, t=0.0, sigma=1.2):
    return np.exp(-(x / sigma) ** 2) * (1.0 + 0.15 * np.sin(2.0 * t))

def frequency_term(t=0.0):
    return 1.0 + 0.2 * np.cos(3.0 * t)

def efcl_potential(x, t=0.0, lam1=0.35, lam2=0.10):
    """Original (unphysical at infinity) version, kept for reference/comparison only."""
    env = env_field(x, t)
    nu = frequency_term(t)
    return lam1 * env + lam2 * (nu ** 2)

def efcl_potential_localized(x, t=0.0, lam1=0.35, lam2=0.10, sigma=1.2):
    """
    Fixed version: the nu^2 term is now modulated by the SAME spatial
    envelope as env_field, so Phi_EFCL -> 0 as x -> +/- infinity.
    This restores a proper asymptotically-free scattering problem.
    """
    envelope = np.exp(-(x / sigma) ** 2)
    env = env_field(x, t, sigma=sigma)
    nu = frequency_term(t)
    return lam1 * env + lam2 * (nu ** 2) * envelope


def analytic_rect_T(E, V0=1.0, a2=2.0, hbar=1.0, m=1.0):
    if E <= 0:
        return 0.0
    if abs(E - V0) < 1e-12:
        return float(1.0 / (1.0 + (m * V0 * a2**2) / (2 * hbar**2)))
    if E < V0:
        kappa = np.sqrt(2 * m * (V0 - E)) / hbar
        T = 1.0 / (1.0 + (V0**2 * np.sinh(kappa * a2)**2) / (4 * E * (V0 - E)))
    else:
        k = np.sqrt(2 * m * (E - V0)) / hbar
        T = 1.0 / (1.0 + (V0**2 * np.sin(k * a2)**2) / (4 * E * (E - V0)))
    return float(T)


def logderiv_T(V_func, E, x_L=-8.0, x_R=8.0, hbar=1.0, m=1.0, **Vkwargs):
    """
    Exact transmission via the Riccati / logarithmic-derivative method
    (unconditionally stable — no exponential blow-up, unlike naive
    amplitude transfer matrices).

    y = psi'/psi obeys  y' = (2m/hbar^2)(V(x)-E) - y^2

    Start at x_R with the pure-outgoing (transmitted) boundary condition
    y(x_R) = i*k_out, integrate backward to x_L, then match to
    incident + reflected plane waves there. Because the potential is
    real, flux conservation gives T = 1 - |r|^2 exactly (no need to
    separately track wavefunction amplitude/normalization).
    """
    V_R = V_func(np.array([x_R]), **Vkwargs)[0]
    V_L = V_func(np.array([x_L]), **Vkwargs)[0]
    k_out = np.sqrt(2 * m * complex(E - V_R)) / hbar
    k_in = np.sqrt(2 * m * complex(E - V_L)) / hbar
    if k_out == 0:
        k_out = 1e-9
    if k_in == 0:
        k_in = 1e-9

    def rhs(x, y_complex):
        y = y_complex[0] + 1j * y_complex[1]
        Vx = V_func(np.array([x]), **Vkwargs)[0]
        dy = (2 * m / hbar**2) * (Vx - E) - y**2
        return [dy.real, dy.imag]

    y0 = [np.real(1j * k_out), np.imag(1j * k_out)]

    sol = solve_ivp(
        rhs, [x_R, x_L], y0, method="RK45",
        rtol=1e-10, atol=1e-12, dense_output=False, max_step=(x_R - x_L) / 4000
    )
    if not sol.success:
        return np.nan, np.nan

    y_L = sol.y[0, -1] + 1j * sol.y[1, -1]

    # reflection coefficient from matching psi = e^{ikx} + r e^{-ikx} at x_L
    r = np.exp(2j * k_in * x_L) * (1j * k_in - y_L) / (1j * k_in + y_L)
    R = np.abs(r) ** 2
    T = 1.0 - R
    return float(T), float(R)


if __name__ == "__main__":
    print("=" * 78)
    print("VALIDATION: log-derivative method vs analytic sharp rectangular barrier")
    print("=" * 78)
    for E in [0.05, 0.3, 0.6, 0.9, 1.0, 1.2, 1.5, 1.7]:
        T_ld, R_ld = logderiv_T(V_barrier_smooth, E, x_L=-8, x_R=8,
                                 V0=1.0, a=1.0, w=0.01)
        T_an = analytic_rect_T(E, V0=1.0, a2=2.0)
        print(f"E={E:.2f} | T_logderiv={T_ld:.8e} | T_analytic={T_an:.8e} | "
              f"rel.err={abs(T_ld-T_an)/max(T_an,1e-30):.3e} | R+T check={R_ld+T_ld:.10f}")

    print()
    print("=" * 78)
    print("Smooth barrier (w=0.12), standard vs EFCL-corrected: log-derivative reference")
    print("=" * 78)
    print(f"{'E':>6} {'T_std_ref':>16} {'T_efcl_ref':>16} {'suppression ratio':>18}")
    for E in [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]:
        Ts, _ = logderiv_T(V_barrier_smooth, E, x_L=-8, x_R=8, V0=1.0, a=1.0, w=0.12)

        def V_efcl_func(x, **kw):
            return V_barrier_smooth(x, V0=1.0, a=1.0, w=0.12) + efcl_potential(x, t=0.0, lam1=0.35, lam2=0.10)

        Te, _ = logderiv_T(V_efcl_func, E, x_L=-8, x_R=8)
        print(f"{E:6.2f} {Ts:16.8e} {Te:16.8e} {Te/max(Ts,1e-300):18.6f}")
