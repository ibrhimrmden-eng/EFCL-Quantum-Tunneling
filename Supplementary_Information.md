# Supplementary Information

## Numerical Scaling Collapse and Systematic Exclusion Analysis for Localized EFCL Corrections in Quantum Tunneling

This document consolidates material referenced but not reproduced in full in the main text: the complete hypothesis-exclusion chain, additional sensitivity checks, and an inventory of diagnostic scripts retained for transparency.

---

## S1. Full hypothesis-exclusion chain

The main text (Discussion, §5.2) summarizes the systematic exclusion analysis. The complete sequence of hypotheses tested, in chronological order, is reproduced here for full transparency.

| # | Hypothesis | Test | Outcome |
|---|---|---|---|
| 1 | Grid error explains the observed suppression | Explicit convergence check (N=500→32000) | Rejected after fix; full convergence confirmed |
| 2 | Perturbative WKB predicts the ln(ratio) curvature | Direct analytic-vs-numeric comparison at E=0.3 (deep tunneling) | Rejected — 36% error in the linear coefficient, 5× error in the quadratic coefficient; dominated by near-turning-point sensitivity |
| 3 | An effective rectangular barrier (V0_eff(λ), L_eff(λ)) explains EFCL behavior generally | Six-parameter fit at one energy (error <2×10⁻⁵), then tested for generalization to a different energy without refitting | Rejected — 30–41% error when transferred from E=0.9 to E=0.5; the excellent local fit was interpolation, not physical structure |
| 4 | The classical action S(E) is more regular than L_eff | Computed ΔS(λ,E)=S_EFCL−S₀ across six energies | Partially rejected — more regular than L_eff, but exhibits a genuine non-analytic threshold near/above the barrier top from quantum over-barrier reflection, which first-order WKB entirely misses |
| 5 | R(E,λ) is separable as F(E)·G(λ) | Coefficient-of-variation test across a full (E,λ) grid | Rejected — CV rises from 2.7% (small λ) to 12.6% (λ=0.40): genuine E–λ coupling |
| 6 | A kink exists in d²R/dλ² | Re-tested with a **uniform** λ grid and proper central differences | Rejected — the apparent kink was a pure artifact of the original non-uniform grid |
| 7 | **Scaling collapse of R(E,λ) onto λ/λ\*(E)** | Direct visual and quantitative collapse check | **Survived** — spread <8% at λ/λ\*=0.5, <1% at λ/λ\*=1.5 |
| 8 | λ\*(E) is a simple quadratic function of E | Fit on E∈[0.2,1.1], tested out-of-sample on three physically distinct regions | Fails catastrophically above the barrier top (RMSE ×85 vs. in-sample), but survives remarkably well through the tunneling→over-barrier crossover region |
| 9 | λ\*(E) reduces to a single semiclassical quantity (E, S₀, L, V̄) | Pearson/Spearman/RMSE comparison | S₀ strongest (r=−0.998) but not established as causal; none fully explanatory |
| 10 | χ=λ·S₀(E) is the natural collapse variable | Direct constancy check of χ\*=λ\*(E)·S₀(E) across E | Rejected — varies by a factor of 8.95× (CV=43.9%) |
| 11 | χ=λ·S₀(E)/L(E) is a better composite variable | Same constancy check | Rejected — improves only marginally (factor 5.7×, CV=36.5%) |
| 12 | λ\*(E) reduces to the scattering phase shift δ(E) or its transforms | Pearson/Spearman/RMSE vs. δ, tan δ, sin δ, cos δ | Rejected — δ(E) is non-monotonic (U-shaped, minimum near E≈0.5–0.6) while λ\*(E) is strictly monotonic, so correlation is structurally weak (Spearman ρ=0.38) |
| 13 | The Wigner time delay dδ/dE explains λ\*(E) | Same correlation battery | Weak — monotonic and better than raw δ (r=0.884), but still weaker than S₀ |
| 14 | The curvature of δ(E) tracks the acceleration of λ\*(E) | d²δ/dE² vs. dλ\*/dE, feature-location comparison | Rejected — the sharp feature in δ'' sits at very low E (deep tunneling); the real acceleration in λ\*(E) occurs near the barrier top. No positional match |
| 15 | The real-axis Jost-function magnitude \|F(k)\| explains λ\*(E) | Correlation test; also compared to T₀(E) directly | Rejected — \|F(k)\|=1/√T₀ by construction, so it carries no information independent of T₀/S₀, and performs worse than S₀ numerically |

**Conclusion carried into the main text:** among all fifteen hypotheses tested, only the scaling collapse itself (hypothesis 7) survives every robustness check applied to it (threshold invariance, barrier-geometry robustness, envelope-width sensitivity). No single first-order semiclassical or single-channel scattering quantity explains λ\*(E).

---

## S2. Additional sensitivity and robustness detail

### S2.1 Numerical solver tolerance (extends Table 7)
Full comparison at three representative energies (E=0.3, 0.7, 0.95), varying `rtol`/`atol`/max-step simultaneously:

| Setting | Max % difference in λ\* vs. baseline |
|---|---|
| Looser (rtol=1e-8, atol=1e-10, 4× fewer steps) | <3×10⁻⁶ % |
| Tighter (rtol=1e-12, atol=1e-14, 2× more steps) | <3×10⁻⁸ % |

The Riccati/log-derivative method shows no meaningful sensitivity to solver tolerance across two orders of magnitude of relaxation — consistent with its theoretical stability property (Section 2.2 of the main text).

### S2.2 Preliminary Gaussian-barrier check (extends §5.3, Limitations)
As a first probe of whether the scaling phenomenon is tied to the specific tanh barrier shape, λ\*(E) was recomputed for a Gaussian barrier V(x)=V₀exp(−x²/2s²) at three energies:

| E | λ\*(E), Gaussian barrier | λ\*(E), tanh barrier (reference) |
|---|---|---|
| 0.2 | 0.2013 | 0.2300 |
| 0.5 | 0.2580 | 0.2694 |
| 0.8 | 0.3572 | 0.3124 |

The absolute values differ (expected, given the different barrier shape), but the qualitative monotonic trend is preserved. This is a preliminary qualitative observation only; it does not re-establish threshold invariance, collapse quality, or separability for the Gaussian family, and should not be read as a validated cross-family result.

### S2.3 Envelope-profile sensitivity (extends Table 6)
Beyond barrier width and height, the coupling profile's own spatial parameters were tested:

| Variation | CV of λ\*(E) ratio | Interpretation |
|---|---|---|
| Envelope width σ (1.0→1.4) | 5.96–6.53% | Sensitive, with a systematic (non-random) trend across energy — physically expected, since σ controls the spatial overlap between the correction and the barrier |
| Envelope center x₀ (−0.3→+0.3) | 1.96% | Robust; ratios at x₀=−0.3 and x₀=+0.3 are numerically identical, as required by the barrier's left–right symmetry about x=0 |

---

## S3. Diagnostic / earlier-stage scripts (not part of the final canonical pipeline)

These scripts are retained in the repository for transparency — they document errors that were found and corrected, or intermediate hypotheses that were tested and excluded (see S1 above) — but are **not** inputs to any table or figure in the main text.

| Script | Role |
|---|---|
| `efcl_test.py` | Original (buggy) model; demonstrates the divergent gradient term later corrected |
| `efcl_v2.py` | First corrected version (smooth barrier); superseded by `logderiv_ref.py` |
| `transfer_matrix_ref.py` | Naive transfer-matrix solver; demonstrated numerically unstable, explicitly not used as the reference method |
| `rigor_checks.py` | Null-limit, λ-sensitivity, and time-dependence sanity checks (hypothesis chain, early stage) |
| `action_deltaS_test.py` | ΔS(λ,E) computation underlying hypothesis 4 in S1 |
| `separability_scaling_test.py` | Separability and initial (non-uniform-grid) collapse test underlying hypotheses 5–7 in S1 |
| `kink_retest_uniform.py` | Uniform-grid re-test that rejected the spurious kink (hypothesis 6 in S1) |
| `verify_effective_barrier_claim.py` | Effective-rectangular-barrier fit and cross-energy generalization test underlying hypothesis 3 in S1 |
| `universal_shift_test.py` | dL_eff/dλ instability check supporting rejection of hypothesis 3 |
| `lambda_star_fit.py` | Quadratic-in-E model fit and out-of-sample RMSE underlying hypothesis 8 in S1 |

Full narrative detail for each is given in `EFCL_scaling_law_results_documentation.md`, included in the repository.

---

## S4. Repository and continuous integration

All code, data, and figures are available at:
https://github.com/ibrhimrmden-eng/EFCL-Quantum-Tunneling

A continuous-integration workflow (`.github/workflows/reproducibility.yml`) independently re-executes the core computational chain (solver validation, WKB comparison, all four parts of Table 4, Figure 4, Figure 8) on a fresh GitHub-hosted Ubuntu runner on every push, as an additional, machine-independent reproducibility guarantee beyond the manual clean-room test described in the repository's README.
