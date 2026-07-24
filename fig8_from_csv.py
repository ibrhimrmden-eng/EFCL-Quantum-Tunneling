import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

rows = []
for fname in ['table4_semiclassical.csv', 'table4_phase_shift.csv', 'table4_jost.csv']:
    with open(f'data_csv/{fname}') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

label_map = {
    'E_raw': 'E (raw)', 'S0': r'$S_0(E)$', 'L': 'L(E)', 'Vbar': r'$\bar V(E)$',
    'delta': r'$\delta(E)$', 'tan_delta': r'$\tan\delta$', 'wigner_delay': r'$d\delta/dE$',
    'F_mag': r'$|F(k)|$', 'ln_F_mag': r'$\ln|F(k)|$'
}
order = ['E_raw','S0','L','Vbar','delta','tan_delta','wigner_delay','F_mag','ln_F_mag']
rows_sorted = sorted(rows, key=lambda r: order.index(r['quantity']))

names = [label_map[r['quantity']] for r in rows_sorted]
pearson = [float(r['pearson_r']) for r in rows_sorted]
rmse = [float(r['rmse']) for r in rows_sorted]

fig, ax = plt.subplots(figsize=(8,5))
y_pos = np.arange(len(names))
colors = ['#2ca02c' if abs(p) > 0.95 else '#ff7f0e' for p in pearson]
ax.barh(y_pos, rmse, color=colors)
ax.set_yticks(y_pos)
ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_xlabel('RMSE of best single-variable linear fit to ' + r'$\lambda^*(E)$')
ax.set_title('Figure 8: Candidate-observable correlation map (from data_csv/table4_*.csv)')
for i,(p,r) in enumerate(zip(pearson,rmse)):
    ax.text(r+0.0008, i, f'r={p:+.3f}', va='center', fontsize=8)
plt.tight_layout()
plt.savefig('figures/Figure8_correlation_map.png', dpi=600)
plt.close()
print("Figure 8 now reads strictly from table4_*.csv files (true canonical pipeline)")
