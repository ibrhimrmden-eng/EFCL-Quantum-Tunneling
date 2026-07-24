import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 11, 'font.family':'serif'})

def read_csv(fname):
    with open(f'data_csv/{fname}') as f:
        return list(csv.DictReader(f))

t4_semi = {r['quantity']: r for r in read_csv('table4_semiclassical.csv')}
t4_phase = {r['quantity']: r for r in read_csv('table4_phase_shift.csv')}
t4_jost = {r['quantity']: r for r in read_csv('table4_jost.csv')}

cv_lamS0 = cv_lamS0L = None
with open('data_csv/table4_composite.csv') as f:
    for row in csv.reader(f):
        if row and row[0] == 'CV_lamS0_percent':
            cv_lamS0 = float(row[1])
        if row and row[0] == 'CV_lamS0_over_L_percent':
            cv_lamS0L = float(row[1])

import numpy as np
E5, r05, r10, r15, r20 = [], [], [], [], []
with open('data_csv/table5_threshold_invariance.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        E5.append(float(row['E']))
        r05.append(float(row['lambda_star_R=-0.5']))
        r10.append(float(row['lambda_star_R=-1.0']))
        r15.append(float(row['lambda_star_R=-1.5']))
        r20.append(float(row['lambda_star_R=-2.0']))
r05,r10,r15,r20 = map(np.array,(r05,r10,r15,r20))
cv_thresh = np.mean([ (r05/r10).std()/(r05/r10).mean(), (r15/r10).std()/(r15/r10).mean(), (r20/r10).std()/(r20/r10).mean() ])*100

def cv_from_csv(fname, ref_col, other_cols):
    rows = read_csv(fname)
    ref = np.array([float(r[ref_col]) for r in rows])
    cvs = []
    for col in other_cols:
        other = np.array([float(r[col]) for r in rows])
        ratio = other/ref
        cvs.append(ratio.std()/ratio.mean())
    return np.mean(cvs)*100

cv_width = cv_from_csv('table6_barrier_width.csv', 'lambda_star_w=0.12', ['lambda_star_w=0.08','lambda_star_w=0.2'])
cv_V0 = cv_from_csv('table6_barrier_V0.csv', 'lambda_star_V0=1.0', ['lambda_star_V0=0.8','lambda_star_V0=1.2'])
cv_sigma = cv_from_csv('table6_envelope_sigma.csv', 'lambda_star_sigma=1.2', ['lambda_star_sigma=1.0','lambda_star_sigma=1.4'])

print(f"cv_thresh={cv_thresh:.2f}%  cv_width={cv_width:.2f}%  cv_V0={cv_V0:.2f}%  cv_sigma={cv_sigma:.2f}%")

categories = [
    ('Semiclassical\n(E, S0, L, Vbar)', max(abs(float(t4_semi['S0']['pearson_r'])), 0), 'Excluded (best: S0)'),
    ('Composite\n(lam*S0, lam*S0/L)', cv_lamS0, f'Excluded (CV {cv_lamS0:.0f}-{cv_lamS0L:.0f}%)'),
    ('Scattering phase\n(delta family)', abs(float(t4_phase['delta']['spearman_rho'])), 'Excluded (non-monotonic)'),
    ('Wigner delay', abs(float(t4_phase['wigner_delay']['pearson_r'])), 'Weak (worse than S0)'),
    ('Jost magnitude\n(real axis)', abs(float(t4_jost['F_mag']['pearson_r'])), 'Excluded (= T0 by definition)'),
    ('Threshold choice', cv_thresh, f'Robust (CV {cv_thresh:.1f}%)'),
    ('Barrier width w', cv_width, f'Robust (CV {cv_width:.2f}%)'),
    ('Barrier height V0', cv_V0, f'Robust w/ rescaling (CV {cv_V0:.1f}%)'),
    ('Envelope sigma', cv_sigma, f'Sensitive (CV {cv_sigma:.1f}%)'),
]

colors9 = ['#d62728','#d62728','#d62728','#ff7f0e','#d62728','#2ca02c','#2ca02c','#2ca02c','#ff7f0e']

fig, ax = plt.subplots(figsize=(9,5.5))
y = np.arange(len(categories))[::-1]
ax.barh(y, [1]*len(categories), color=colors9, height=0.6)
ax.set_yticks(y)
ax.set_yticklabels([c[0] for c in categories], fontsize=9)
for yi, c in zip(y, categories):
    ax.text(0.5, yi, c[2], ha='center', va='center', fontsize=8, color='white', fontweight='bold')
ax.set_xticks([])
ax.set_title('Figure 9: Systematic exclusion map (from data_csv/table4-6*.csv)')
legend_elems = [mpatches.Patch(color='#d62728', label='Excluded'),
                mpatches.Patch(color='#ff7f0e', label='Weak / partial'),
                mpatches.Patch(color='#2ca02c', label='Robust (survives test)')]
ax.legend(handles=legend_elems, loc='lower right', frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig('figures/Figure9_exclusion_map.png', dpi=600)
plt.close()
print("Figure 9 now reads strictly from CSV files; verdict labels are the only remaining editorial content")
