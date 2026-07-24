import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

# ---------------- Figure 8: Correlation map ----------------
candidates = ['E (raw)', r'$S_0(E)$', 'L(E)', r'$\bar V(E)$', r'$\lambda S_0$', r'$\lambda S_0/L$',
              r'$\delta(E)$', r'$\tan\delta$', r'$d\delta/dE$', r'$|F(k)|$', r'$\ln|F(k)|$', r'$T_0(E)$']
pearson = [0.998, -0.998, -0.978, 0.908, np.nan, np.nan, -0.516, -0.551, 0.884, -0.731, -0.938, 0.971]
rmse    = [0.00334, 0.00363, 0.01109, 0.02202, np.nan, np.nan, 0.0451, 0.0439, 0.0246, 0.0359, 0.0183, 0.0126]
cv_percent = [None,None,None,None,43.9,36.5,None,None,None,None,None,None]

fig, ax = plt.subplots(figsize=(8,5))
y_pos = np.arange(len(candidates))
rmse_plot = [r if not np.isnan(r) else 0 for r in rmse]
def pick_color(p):
    if np.isnan(p):
        return '#999999'
    return '#2ca02c' if abs(p) > 0.95 else '#ff7f0e'
colors = [pick_color(p) for p in pearson]
bars = ax.barh(y_pos, rmse_plot, color=colors)
ax.set_yticks(y_pos)
ax.set_yticklabels(candidates)
ax.invert_yaxis()
ax.set_xlabel('RMSE of best single-variable linear fit to ' + r'$\lambda^*(E)$')
ax.set_title(r'Figure 8: Candidate-observable correlation map (lower RMSE = stronger)')
for i,(p,r) in enumerate(zip(pearson,rmse)):
    if not np.isnan(p):
        ax.text(r+0.0008, i, f'r={p:+.3f}', va='center', fontsize=8)
    else:
        ax.text(0.0008, i, 'not constant (CV=%.1f%%)'%cv_percent[i], va='center', fontsize=8, color='#d62728')
plt.tight_layout()
plt.savefig('figures/Figure8_correlation_map.png', dpi=600)
plt.close()
print("Figure 8 done")

# ---------------- Figure 9: Systematic exclusion map ----------------
categories = ['Semiclassical\n(E, S0, L, Vbar)', 'Composite\n(lam*S0, lam*S0/L)', 'Scattering phase\n(delta, tan/sin/cos)',
              'Wigner delay\n(d(delta)/dE)', 'Jost magnitude\n(real axis |F(k)|)', 'Threshold choice\n(R=-0.5..-2.0)',
              'Barrier width w\n(0.08-0.20)', 'Barrier height V0\n(0.8-1.2)', 'Envelope sigma\n(1.0-1.4)']
status = ['Excluded (best: S0, r=-0.998,\nnot causal)', 'Excluded (CV 36-44%)', 'Excluded (Spearman=0.38)',
          'Weak (r=0.88, worse than S0)', 'Excluded (no new info,\n= T0 by definition)', 'Robust (CV 2-3.5%)',
          'Robust (CV<1%)', 'Robust w/ rescaling\n(CV<2%)', 'Sensitive (CV~6%,\nsystematic trend)']
colors9 = ['#d62728','#d62728','#d62728','#ff7f0e','#d62728','#2ca02c','#2ca02c','#2ca02c','#ff7f0e']

fig, ax = plt.subplots(figsize=(9,5.5))
y = np.arange(len(categories))[::-1]
ax.barh(y, [1]*len(categories), color=colors9, height=0.6)
ax.set_yticks(y)
ax.set_yticklabels(categories, fontsize=9)
for yi, s in zip(y, status):
    ax.text(0.5, yi, s, ha='center', va='center', fontsize=8, color='white', fontweight='bold')
ax.set_xticks([])
ax.set_title('Figure 9: Systematic exclusion map for ' + r'$\lambda^*(E)$' + ' candidate explanations')
legend_elems = [mpatches.Patch(color='#d62728', label='Excluded'),
                mpatches.Patch(color='#ff7f0e', label='Weak / partial'),
                mpatches.Patch(color='#2ca02c', label='Robust (survives test)')]
ax.legend(handles=legend_elems, loc='lower right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig('figures/Figure9_exclusion_map.png', dpi=600)
plt.close()
print("Figure 9 done")

# ---------------- Figure 10: Workflow diagram ----------------
steps = ["Original model\n(sharp barrier +\ndivergent grad. term)",
         "Error detected\n(1/dx divergence;\nnon-localized coupling)",
         "Correction\n(smooth barrier +\nlocalized envelope)",
         "Validation\n(log-derivative solver,\n<0.03% vs analytic)",
         "Scaling collapse\nR(E,\u03bb) \u2192 H(\u03bb/\u03bb*(E))",
         "Systematic exclusion\n(11 candidates tested,\nnone fully explanatory)"]

fig, ax = plt.subplots(figsize=(10,3.2))
n = len(steps)
box_w, box_h = 1.55, 1.3
xs = np.linspace(0, (n-1)*1.9, n)
for i,(x0,txt) in enumerate(zip(xs, steps)):
    color = '#e8f0fe' if i not in (1,) else '#fde8e8'
    box = FancyBboxPatch((x0-box_w/2, -box_h/2), box_w, box_h,
                          boxstyle="round,pad=0.08", linewidth=1.3,
                          edgecolor='#333333', facecolor=color)
    ax.add_patch(box)
    ax.text(x0, 0, txt, ha='center', va='center', fontsize=7.8)
    if i < n-1:
        arrow = FancyArrowPatch((x0+box_w/2, 0), (xs[i+1]-box_w/2, 0),
                                 arrowstyle='-|>', mutation_scale=14, color='#333333', lw=1.3)
        ax.add_patch(arrow)

ax.set_xlim(xs[0]-box_w, xs[-1]+box_w)
ax.set_ylim(-1.2, 1.2)
ax.axis('off')
ax.set_title('Figure 10: Study workflow', fontsize=12)
plt.tight_layout()
plt.savefig('figures/Figure10_workflow.png', dpi=600)
plt.close()
print("Figure 10 done")
