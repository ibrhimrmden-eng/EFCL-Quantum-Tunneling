Dear Editor,

We are pleased to submit our manuscript, "Numerical Scaling Collapse and Systematic Exclusion Analysis for Localized EFCL Corrections in Quantum Tunneling," for consideration as a Paper in Journal of Physics A: Mathematical and Theoretical.

**The problem.** Proposals to modify quantum-tunneling transmission through an additive, environment- or frequency-dependent correction to a potential barrier require a numerically reliable reference solver, since the corrections under study are frequently smaller than the errors introduced by an inadequate numerical method. In the course of testing one such correction (the Environment-Frequency Confinement Law, EFCL), we identified and corrected two methodological defects in its original numerical formulation — a spatial derivative term that diverges under grid refinement, and a coupling term that violates the asymptotic-flatness condition required to define a scattering transmission coefficient — and constructed a validated logarithmic-derivative (Riccati) solver, agreeing with the exact analytic transmission coefficient for a rectangular barrier to better than 0.03%.

**What is new.** Using this reference, we show that the standard WKB approximation, including its leading-order perturbative expansion in the coupling strength, fails systematically to reproduce the corrected quantity, with discrepancies of 30–90% near and above the barrier top due to quantum over-barrier reflection. More significantly, we report a robust numerical scaling collapse: the log-suppression ratio R(E,λ) collapses approximately onto a single-variable function of λ/λ*(E), where λ*(E) is a critical coupling scale. This collapse is stable under substantial variation of the barrier's geometric parameters and is approximately independent of the threshold used to define λ*(E). We then conduct a systematic exclusion analysis, ruling out eleven candidate first-order semiclassical and single-channel scattering observables (classical action, turning-point separation, average barrier height, scattering phase shift and its transformations, Wigner delay, and Jost-function magnitude) as adequate explanations for λ*(E).

**Why this fits Journal of Physics A.** The manuscript's central contribution is a numerically exact and independently reproducible finding — a nontrivial scaling relation whose physical origin resists explanation by any tested first-order semiclassical quantity — presented alongside a complete, versioned, and continuously-tested code repository. We believe this combination of a genuine mathematical-physics observation with an unusually thorough negative-result (exclusion) analysis is well suited to the journal's scope and readership, and is presented explicitly as a numerical discovery accompanied by a systematic exclusion map, rather than as a claim of a fully derived physical theory.

**Reproducibility.** All code, data, and figures are archived in a public repository with a continuous-integration workflow that independently re-executes the core computational pipeline on every update (https://github.com/ibrhimrmden-eng/EFCL-Quantum-Tunneling).

This manuscript has not been published previously and is not under consideration elsewhere. We have no conflicts of interest to declare. We believe it will be of interest to researchers working on semiclassical approximations, quantum scattering theory, and numerical methods in quantum mechanics.

Thank you for your consideration.

Sincerely,
Ibrahim Ramadan Al-Shtaiwi
Independent Researcher
