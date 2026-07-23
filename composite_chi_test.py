"""
CANONICAL SCRIPT for the composite chi=lambda*S0 and chi=lambda*S0/L constancy
checks (Table 4's composite rows). Re-derives lambda*(E), S0(E), L(E) directly
rather than reusing values transcribed elsewhere.
"""
import numpy as np
import csv
from lambda_star_extrapolation_test import find_lambda_star
from physical_correlation_test import S0_L_Vbar

E_list = [0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99]

rows = []
for E in E_list:
    ls, T0 = find_lambda_star(E)
    S0, L, Vbar = S0_L_Vbar(E)
    chi_S0 = ls * S0
    chi_S0_L = ls * S0 / L
    rows.append([E, ls, S0, L, chi_S0, chi_S0_L])
    print(f"E={E:.2f}  lambda*={ls:.5f}  S0={S0:.5f}  L={L:.5f}  chi(lam*S0)={chi_S0:.5f}  chi(lam*S0/L)={chi_S0_L:.5f}")

chi_S0_arr = np.array([r[4] for r in rows])
chi_S0_L_arr = np.array([r[5] for r in rows])
cv1 = chi_S0_arr.std()/chi_S0_arr.mean()
cv2 = chi_S0_L_arr.std()/chi_S0_L_arr.mean()
print(f"\nCV(lambda*S0) = {cv1*100:.1f}%")
print(f"CV(lambda*S0/L) = {cv2*100:.1f}%")

with open('data_csv/table4_composite.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['E','lambda_star','S0','L','chi_lamS0','chi_lamS0_over_L'])
    writer.writerows(rows)
    writer.writerow([])
    writer.writerow(['CV_lamS0_percent', cv1*100])
    writer.writerow(['CV_lamS0_over_L_percent', cv2*100])
print("Saved: data_csv/table4_composite.csv")
