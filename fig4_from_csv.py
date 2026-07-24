import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

data = {}
with open('data_csv/figure4_scaling_collapse.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        E = float(row['E'])
        data.setdefault(E, {'lam':[], 'R':[], 'lam_star':None, 'x':[]})
        data[E]['lam'].append(float(row['lambda']))
        data[E]['R'].append(float(row['R']))
        data[E]['lam_star'] = float(row['lambda_star'])
        x = row['lambda_over_lambda_star']
        data[E]['x'].append(float(x) if x != '' else np.nan)

E_list = sorted(data.keys())
colors = plt.cm.viridis(np.linspace(0,1,len(E_list)))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4.5))
for E, c in zip(E_list, colors):
    d = data[E]
    ax1.plot(d['lam'], d['R'], 'o-', color=c, label=f'E={E}', ms=4)
    ax2.plot(d['x'], d['R'], 'o-', color=c, label=f'E={E}', ms=4)

ax1.set_xlabel(r'$\lambda$')
ax1.set_ylabel(r'$R(E,\lambda)=\ln(T_{EFCL}/T_0)$')
ax1.set_title('(a) Raw suppression curves')
ax1.axhline(-1, color='gray', lw=0.7, ls=':')
ax1.legend(frameon=False, fontsize=8)

ax2.set_xlabel(r'$\lambda/\lambda^*(E)$')
ax2.set_ylabel(r'$R(E,\lambda)$')
ax2.set_title(r'(b) Rescaled by $\lambda^*(E)$: scaling collapse')
ax2.axhline(-1, color='gray', lw=0.7, ls=':')
ax2.axvline(1, color='gray', lw=0.7, ls=':')
ax2.legend(frameon=False, fontsize=8)

fig.suptitle(r'Figure 4: Scaling collapse of $R(E,\lambda)$ onto $\lambda/\lambda^*(E)$')
plt.tight_layout()
plt.savefig('figures/Figure4_scaling_collapse.png', dpi=600)
plt.close()
print("Figure 4 regenerated directly from figure4_scaling_collapse.csv")
