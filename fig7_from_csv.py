import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

E_s, lam_s = [], {1.0:[], 1.2:[], 1.4:[]}
with open('data_csv/table6_envelope_sigma.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        E_s.append(float(row['E']))
        lam_s[1.0].append(float(row['lambda_star_sigma=1.0']))
        lam_s[1.2].append(float(row['lambda_star_sigma=1.2']))
        lam_s[1.4].append(float(row['lambda_star_sigma=1.4']))

fig, ax = plt.subplots(figsize=(6,4.2))
for sig, c in zip([1.0,1.2,1.4], ['#1f77b4','#2ca02c','#d62728']):
    ax.plot(E_s, lam_s[sig], 'o-', color=c, label=fr'$\sigma$={sig}')
ax.set_xlabel('Energy E'); ax.set_ylabel(r'$\lambda^*(E)$')
ax.set_title(r'Figure 7: Sensitivity to EFCL envelope width $\sigma$')
ax.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig('figures/Figure7_sigma_sensitivity.png', dpi=600)
plt.close()
print("Figure 7 now reads strictly from CSV (true canonical pipeline)")
