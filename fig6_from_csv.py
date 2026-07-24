import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

E_w, lam_w = [], {0.08:[], 0.12:[], 0.20:[]}
with open('data_csv/table6_barrier_width.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        E_w.append(float(row['E']))
        lam_w[0.08].append(float(row['lambda_star_w=0.08']))
        lam_w[0.12].append(float(row['lambda_star_w=0.12']))
        lam_w[0.20].append(float(row['lambda_star_w=0.2']))

frac_V0, lam_V0 = [], {0.8:[], 1.0:[], 1.2:[]}
with open('data_csv/table6_barrier_V0.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        frac_V0.append(float(row['frac_E_over_V0']))
        lam_V0[0.8].append(float(row['lambda_star_V0=0.8']))
        lam_V0[1.0].append(float(row['lambda_star_V0=1.0']))
        lam_V0[1.2].append(float(row['lambda_star_V0=1.2']))

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4.2))
for wv, c in zip([0.08,0.12,0.20], ['#1f77b4','#2ca02c','#d62728']):
    ax1.plot(E_w, lam_w[wv], 'o-', color=c, label=f'w={wv}')
ax1.set_xlabel('Energy E'); ax1.set_ylabel(r'$\lambda^*(E)$')
ax1.set_title('(a) Edge-sharpness variation')
ax1.legend(frameon=False, fontsize=9)

for V0v, c in zip([0.8,1.0,1.2], ['#1f77b4','#2ca02c','#d62728']):
    ax2.plot(frac_V0, lam_V0[V0v], 's-', color=c, label=f'V0={V0v}')
ax2.set_xlabel('E / V0'); ax2.set_ylabel(r'$\lambda^*(E)$')
ax2.set_title('(b) Barrier-height variation')
ax2.legend(frameon=False, fontsize=9)

fig.suptitle('Figure 6: Robustness of the collapse to barrier geometry')
plt.tight_layout()
plt.savefig('figures/Figure6_barrier_robustness.png', dpi=600)
plt.close()
print("Figure 6 now reads strictly from CSV (true canonical pipeline)")
