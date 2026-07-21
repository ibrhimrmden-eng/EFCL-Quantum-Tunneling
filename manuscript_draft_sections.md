# Manuscript Draft Sections
### Numerical Scaling Collapse and Systematic Exclusion Analysis for Localized EFCL Corrections in Quantum Tunneling

---

## Abstract

We present a numerically exact investigation of a proposed environment-dependent correction, Φ_EFCL, to one-dimensional quantum tunneling through a smooth potential barrier. Two methodological defects in the original formulation are identified and corrected: a spatial derivative term that diverges under grid refinement, and a spatially non-decaying coupling term that violates the asymptotic-flatness condition required to define a scattering transmission coefficient. All numerical results are computed against a validated reference — a logarithmic-derivative (Riccati) solution of the time-independent Schrödinger equation, agreeing with the exact analytic transmission coefficient for a rectangular barrier to better than 0.03% and satisfying flux conservation (R+T=1) to at least ten decimal digits.

Using this reference, we show that the standard Wentzel–Kramers–Brillouin (WKB) approximation, including its leading-order perturbative expansion in the coupling strength λ, fails systematically to reproduce the corrected quantity, with discrepancies of 30–90% near and above the barrier top due to quantum over-barrier reflection, an effect absent from the primitive WKB exponent. We report a robust numerical scaling collapse: the log-suppression ratio R(E,λ) = ln[T_EFCL(E,λ)/T₀(E)] collapses approximately onto a single-variable function of λ/λ*(E), where λ*(E) is defined by R(E,λ*)=−1. The collapse is approximately independent of this threshold choice (systematic deviation of 2–3.5% across four tested thresholds) and stable under substantial variation of the barrier's geometric parameters (width and height, CV<2%), with a preliminary indication that it is not restricted to the specific barrier shape studied. A systematic exclusion analysis shows that λ*(E) is not determined by any single first-order semiclassical or single-channel scattering observable tested — including the classical action, turning-point separation, average barrier height, the scattering phase shift, its derivatives, or the real-axis Jost-function magnitude — although the classical action S₀(E) is the strongest tested predictor. We frame these results as a numerical discovery and a systematic exclusion map, explicitly distinct from a first-principles derivation of the observed scaling.

---

## 1. Introduction

Proposals to modify quantum-tunneling transmission through an additive, environment- or frequency-dependent correction to the potential barrier are a natural extension of density-dependent modifications explored in other contexts, including the Environment-Frequency Confinement Law (EFCL) framework under investigation here [Al-Shtaiwi, 2026, Zenodo DOI:10.5281/zenodo.19936345]. Testing such a correction numerically requires (i) a barrier and coupling profile that are well-posed as a scattering problem, and (ii) a transmission solver whose accuracy and stability are independently established, since the corrections under study are frequently smaller than the errors introduced by an inadequate numerical method.

The WKB approximation is the standard starting point for tunneling calculations, but its limitations are well documented. The primitive (non-uniform) WKB exponent is known to break down in the vicinity of classical turning points, where a uniform treatment based on Airy-function connection formulas is required [Berry & Mount, Rep. Prog. Phys. 35, 315 (1972)]. Independently, WKB in its simplest form does not capture quantum reflection above a barrier — a purely quantum-mechanical effect in which a particle with energy exceeding the barrier maximum is reflected with non-negligible probability, absent from the classical-turning-point picture on which the primitive WKB exponent is built [Jaffe, arXiv:0910.4564]. Both limitations are directly relevant to any perturbative or semiclassical treatment of an energy-dependent correction whose associated coupling strength can, for some energies, push the effective barrier above or below the incident energy.

For a numerically reliable reference beyond WKB, we adopt the logarithmic-derivative (Riccati) method, a well-established and numerically stable alternative to direct propagation of transfer matrices, which are known to suffer catastrophic loss of precision for optically thick barriers due to the multiplicative accumulation of exponentially growing and decaying solution components [Johnson, J. Chem. Phys. 67, 4086 (1977)]. Because the log-derivative method propagates a bounded nonlinear (Riccati) variable rather than the wavefunction amplitude itself, it remains stable across the full range of barrier heights and widths considered in this work.

Beyond the semiclassical and single-channel scattering observables examined here (energy, classical action, turning-point separation, average potential, scattering phase shift and its derivatives, and the real-axis magnitude of the Jost function), a more complete scattering-theoretic description requires the analytic continuation of the Jost function to complex momentum, including the identification of its zeros and associated resonance poles [Newton, *Scattering Theory of Waves and Particles*, 2nd ed., Springer-Verlag, 1982]. This complex-plane analysis constitutes a distinct research program and is explicitly deferred to future work (Section 5).

**Contributions of this work:**
1. Identification and correction of two methodological errors in the original EFCL-tunneling formulation (a divergent grid-dependent term and a non-localized coupling term violating asymptotic flatness).
2. Construction and validation of a numerically exact log-derivative reference solver, used as the sole basis for all subsequent quantitative claims.
3. A quantitative demonstration that WKB — both in its primitive form and in a leading-order perturbative expansion in the coupling strength — fails to reproduce the exact suppression behavior, with the largest discrepancies located at and above the barrier top.
4. Discovery and characterization of a numerical scaling collapse in the coupling-suppression relationship, together with a threshold-invariance check, a barrier-family robustness check, and a numerical-uncertainty (solver-tolerance) check.
5. A systematic exclusion analysis ruling out eleven candidate first-order semiclassical and single-channel scattering observables as explanatory variables for the extracted critical scale λ*(E).

---

## 2. Mathematical Formulation

### 2.1 Barrier and coupling potential
We consider one-dimensional scattering (ħ=m=1) off a smooth barrier

  V_std(x) = (V₀/2) [tanh((x+a)/w) − tanh((x−a)/w)],

parameterized by height V₀, half-width a, and edge-sharpness w. The proposed EFCL correction is

  Φ_EFCL(x,t) = λ₁·Env(x,t) + λ₂·ν(t)² · envelope(x),

  Env(x,t) = exp[−(x/σ)²]·(1+0.15 sin 2t),  ν(t) = 1+0.2 cos 3t,  envelope(x) = exp[−(x/σ)²],

with σ the spatial width of the coupling profile. The second term was reformulated relative to its original definition, which lacked the envelope(x) factor and consequently failed to vanish as x→±∞; this violates the asymptotic-flatness condition required for a well-defined scattering transmission coefficient and produces unphysical (negative) transmission at low energy in the uncorrected form. The effective potential used throughout is

  V_eff(x;λ,t) = V_std(x) + Φ_EFCL(x,t; λ),

with λ denoting the coupling strength (λ₁ or λ₂ as appropriate to the test in question) and t=0 fixed unless otherwise stated.

### 2.2 Reference solver: logarithmic-derivative (Riccati) method
Define y(x) = ψ'(x)/ψ(x). Substitution into the time-independent Schrödinger equation −ψ''+V(x)ψ=Eψ gives the Riccati equation

  y'(x) = 2[V(x)−E] − y(x)².

The solution is initialized at x=x_R (asymptotically flat) with the pure-outgoing boundary condition y(x_R)=ik_out, k_out=√[2(E−V(x_R))], integrated backward to x=x_L, and matched to an incident-plus-reflected plane wave ψ=e^{ikx}+re^{−ikx} to extract the reflection amplitude r. Because V(x) is real, flux conservation is exact: R+T=1 with R=|r|², independent of any amplitude normalization. This method avoids the numerical instability of amplitude-based transfer-matrix propagation, since y(x) remains bounded rather than tracking exponentially growing and decaying wavefunction components separately. Numerical integration uses adaptive Runge–Kutta (RK45) with rtol=10⁻¹⁰, atol=10⁻¹². Sensitivity to these tolerances is negligible: relaxing them by two orders of magnitude changes extracted quantities by less than 10⁻⁵%.

Validation against the exact analytic transmission coefficient for a sharp rectangular barrier,

  T(E) = [1 + V₀² sinh²(κL)/(4E(V₀−E))]⁻¹  (E<V₀),  κ=√[2(V₀−E)],

gives relative agreement better than 3×10⁻⁴ across the tested energy range, with flux conservation R+T=1 satisfied to at least 10 decimal digits.

### 2.3 Suppression ratio and critical scale
For a given energy E and coupling strength λ, define

  R(E,λ) = ln[ T_EFCL(E,λ) / T₀(E) ],  T₀(E) ≡ T(E; λ=0).

The critical scale λ*(E) is defined implicitly by R(E, λ*(E)) = −1, i.e., the coupling strength at which the corrected transmission is suppressed by a factor of e relative to the bare barrier. Robustness of λ*(E) to this threshold choice is assessed by repeating the extraction at R=−0.5, −1.5, −2.0.

### 2.4 Candidate explanatory quantities
For the deep-tunneling regime (E<V₀, where classical turning points x₁,₂(E) exist for V_std), the following quantities are computed for comparison against λ*(E):

  S₀(E) = ∫_{x₁}^{x₂} √[2(V_std(x)−E)] dx   (classical action),
  L(E) = x₂(E) − x₁(E)          (turning-point separation),
  V̄(E) = (1/L)∫_{x₁}^{x₂} V_std(x) dx    (average barrier height),
  δ(E) = ½ arg[r(E)]           (scattering phase shift, bare barrier),
  dδ/dE                  (Wigner time delay, up to ħ),
  |F(k)| ≡ 1/√T₀(E)            (real-axis Jost-function magnitude, simplest normalization).

Correlation with λ*(E) is assessed via Pearson and Spearman coefficients and the root-mean-square error (RMSE) of the best single-variable linear fit; out-of-sample generalization is assessed by fitting on one energy range and evaluating prediction error on physically distinct, non-overlapping energy regions (deep tunneling, barrier-top crossover, over-barrier).

---

## 5. Discussion and Conclusion

### 5.1 Summary of established results
Within the tested parameter range, three findings are established with high confidence:

1. **The corrected numerical model is well-posed.** After removing the divergent finite-difference term and localizing the coupling potential, the transmission coefficient converges under grid refinement, satisfies the λ=0 null-limit check exactly, and conserves flux to machine precision throughout.
2. **WKB is not a reliable quantitative tool for this problem**, even in its perturbative form and even at energies conventionally considered "safe" for the deep-tunneling regime. The leading-order perturbative expansion of the WKB exponent in λ disagrees with the exact suppression ratio by 30–40% already at E=0.3 (well below the barrier height), and the discrepancy grows to order-of-magnitude near and above the barrier top, where the primitive WKB exponent misses quantum over-barrier reflection entirely.
3. **A robust numerical scaling collapse exists.** R(E,λ)/behavior organizes approximately as a function of the single variable λ/λ*(E), with the collapse quality insensitive to the threshold used to define λ*(E) (2–3.5% systematic deviation) and to the barrier's geometric parameters (width and height, CV<2%), with weaker but non-trivial sensitivity to the spatial width of the coupling profile itself (CV≈6%, physically expected given that this width controls the spatial overlap between the correction and the barrier).

### 5.2 What remains unresolved
The physical origin of λ*(E) is not established. A systematic exclusion analysis — the most extensive component of this work — rules out eleven candidate explanatory variables drawn from first-order semiclassical theory (raw energy, classical action, turning-point separation, average potential, and simple combinations thereof) and single-channel scattering theory (the scattering phase shift and its low-order transformations, the Wigner time delay, and the real-axis Jost-function magnitude) as adequate single-variable predictors of λ*(E). The classical action S₀(E) is the strongest of these (Pearson r=−0.998, RMSE=0.0036), but out-of-sample and cross-quantity tests show this correlation does not extend to a closed-form or causal relationship: a quadratic fit of λ*(E) that performs excellently within its calibration range fails by nearly two orders of magnitude when extrapolated to energies significantly above the barrier top.

We emphasize a distinction that is easily blurred in this kind of analysis: the failure to identify an explanatory variable is a statement about the tested observables, not a statement about the underlying physical model. Nothing in the present results indicates a defect in the EFCL correction as defined (Section 2.1); the exclusion results concern only the adequacy of first-order semiclassical and single-channel descriptors for organizing its energy dependence.

### 5.3 Limitations
(i) All quantitative scaling results derive from a single coupling profile family (Gaussian-enveloped) evaluated primarily on a tanh-shaped barrier; a preliminary check on a Gaussian barrier reproduces the qualitative (monotonic) trend in λ*(E) but does not re-establish threshold invariance or collapse quality for that family. (ii) The exclusion analysis is confined to the deep-tunneling regime (E<V₀) for the semiclassical quantities, since turning points are undefined for the bare barrier above V₀. (iii) The present work does not attempt an analytic derivation of the scaling collapse; it establishes the collapse as a numerical fact and characterizes its robustness and its resistance to simple explanation.

### 5.4 Future work
The natural next step, motivated directly by the exclusion results, is analytic continuation of the scattering solution to complex momentum k, enabling study of the Jost function's zeros and associated resonance poles — a structure not accessible from real-axis observables alone [Newton, 1982]. This is treated as an independent research program rather than an extension of the present numerical study, given the substantially different mathematical machinery required (complex-k propagation, pole-finding, and resonance identification) relative to the real-axis correlational analysis reported here.

### 5.5 Conclusion
We have presented a numerically validated, self-consistent investigation of a localized EFCL correction to one-dimensional quantum tunneling, correcting two methodological errors in the original formulation, establishing the quantitative inadequacy of WKB for this problem, and reporting a robust numerical scaling collapse whose physical origin remains an open question after a systematic exclusion of the most natural first-order candidate explanations. We present this as a numerical discovery accompanied by a systematic exclusion map, rather than as a first-principles theoretical derivation.

---

*Note: This draft intentionally omits a "Related Work" subsection with a broader literature survey beyond the specific technical citations above, since a comprehensive survey of prior EFCL-adjacent or tunneling-correction literature was not part of the verified computational work in this session and should not be fabricated. The author should expand Section 1 with a survey of directly comparable prior proposals (if any exist) before submission.*
