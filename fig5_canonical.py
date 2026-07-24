import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

E_vals = []
data = {t: [] for t in [-0.5,-1.0,-1.5,-2.0]}
with open('data_csv/table5_threshold_invariance.csv') as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        E_vals.append(float(row[0]))
        for i,t in enumerate([-0.5,-1.0,-1.5,-2.0]):
            data[t].append(float(row[i+1]))

E_vals = np.array(E_vals)
ref = np.array(data[-1.0])

fig, ax = plt.subplots(figsize=(6,4.2))
for t in [-0.5,-1.0,-1.5,-2.0]:
    vals = np.array(data[t])
    ax.plot(E_vals, vals/ref, 'o-', label=f'R={t}')
ax.axhline(1.0, color='gray', lw=1, ls=':')
ax.set_xlabel('Energy E')
ax.set_ylabel(r'$\lambda^*(R)/\lambda^*(R=-1)$')
ax.set_title('Figure 5: Threshold invariance of the collapse')
ax.legend(frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig('figures/Figure5_threshold_invariance.png', dpi=600)
plt.close()
print("Figure 5 regenerated directly from table5_threshold_invariance.csv -- exact match with Table 5")
