# Release validation contract

Dated 8 October 2026. These rules were specified after the first package failed a cross-platform replay and before testing this repaired candidate. They are post-hoc engineering rules, **not** a new preregistration or a replacement for the original scientific tolerances.

## Three separate outcomes

1. Scientific reproduction uses the frozen paper claim criteria, with dated amendments and final adjudications. C6 has a narrow strict-tolerance miss; C9 is a contextual cross-dataset statement, not a failed within-DESI claim. A replay PASS does not change either fact.
2. Numerical replay compares deterministic values and stochastic posterior summaries using the contracts below. The original reference outputs and their historical FAIL labels remain intact.
3. Distribution integrity verifies full SHA-256 hashes of every distributed file, excluding the self-referential manifest. Runtime work is confined to `_runs/`. Downloaded data must pass their separately pinned hashes. This verifies consistency, not an external signature or scientific correctness.

## Acceptance criteria

- Deterministic output: text structure and verdict labels must match, with numeric absolute and relative tolerance 1e-9. This parser is restricted to deterministic outputs; it is never applied to chain hashes or sampler diagnostics.
- C9: only displayed fitted nuisance parameters have looser, units-aware tolerances (Ωm 0.0002; H0rd 0.002 km/s). χ² and p-value displays retain the deterministic tolerance.
- C4: deterministic prior comparisons and mock results must replay; each sampler summary must agree with the refined grid within 10% of its smaller 68% interval offset and meet the original paper tolerances. Each parameter must have finite rank-normalised split R-hat < 1.01 and retained steps / autocorrelation time ≥ 50. If the initial 10,000 retained steps fail convergence, repeat once with 20,000 and the same seeds and burn-in. Grid refinement max change ≤ 0.001 and quadrature relative error ≤ 1e-8.
- C8: validate both input identities and both bin selections. For each of Ωm, w0, wa and H0rd, compare each of the three quantiles to the reference within 10% of its smaller reference 68% interval offset, with an additional absolute cap 0.02 for w0. This is a pragmatic cross-platform stability bound, not a Monte Carlo standard-error estimate. Each parameter must meet R-hat < 1.01 and retained steps / autocorrelation time ≥ 50, with quadrature error ≤ 1e-8. Repeat once at 40,000 retained steps if 20,000 fail convergence, preserving seeds and burn-in; record both diagnostics. Failure after the extension fails the run. The historical DR1 rounded-redshift paper criteria are checked separately.
- R-hat uses one fixed walker from each of four independently started ensembles because walkers within an ensemble are coupled. This is a diagnostic with limited coverage, not proof that all posterior modes were explored. Priors and weak wa constraints remain part of the interpretation.
- The suite is complete only with `python3 reproduce.py --all`. CI runs the full suite and negative regression tests. No CI step publishes, creates a release, or records owner approval.

## Interpretation and provenance

Original scientific reference outputs are preserved. Newly calculated posterior summaries, attempt diagnostics and chains are saved separately in the run directory. Chain hashes are provenance identifiers, never numeric quantities to compare across environments. Original expected outcomes in PROTOCOL.md are historical; in particular, its C1 monotonicity expectation and C9 adverse interpretation were corrected during review.

DR1/DR2 comparisons use the same redshift-precision convention. N1 has its own output comparing both priors and the Planck significance; N2 points to the full-precision significance calculation; N4 includes the defining AP identity. DES release tag 1.3 is the archived original 2024 analysis, even though the GitHub tag was published in 2025; it is not the later Dovekie reanalysis.

Approval is a human action on an exact package hash after verification. This builder always emits an unapproved candidate. The current verification session reused earlier review context; it cannot be described as a fresh-conversation clean-room review.
