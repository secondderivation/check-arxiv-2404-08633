> **Historical preregistration and addenda, preserved for audit.** This document records expectations, including expectations later shown wrong. It is not the final verdict. The original C1 “decreasing” expectation was wrong: the ratio increases with Ωm. The original C9 criterion asked the wrong within-dataset question for the paper’s cross-dataset statement; the adverse interpretation was withdrawn. C6 retains its narrow strict-tolerance failure. Final interpretations are in README.md; separate, post-hoc release replay criteria are in VALIDATION.md. References to paper.json, B-BRIEF.md, private reports, and internal workflow stages describe the historical review and are not additional public package inputs. Pinned public data are in datasets.json.

# PROTOCOL: arXiv:2404.08633v3 (Ó Colgáin et al., "Does DESI 2024 Confirm ΛCDM?")

Committed 2026-10-07 before any computation (G2). Changes go in a dated addendum below, also before the run they affect. Anything not planned here is labelled `post_hoc` in the report.

Inputs: `paper.json` (V0 freeze and hashes), `claims/C1–C9-r1.json`, `conventions.json`, `paper_values.json` (printed numbers, transcribed), `data/` (DESI DR1 and DR2 BAO, hashes in `data/SHA256SUMS.txt`).

Roles: this session is **A** (Claude). Everything it produces is **V2 only** (G19). V4 needs a blind **B** run by Codex that commits its own numbers before reading A's outputs (brief: `B-BRIEF.md`, written at V3).

Honesty note: while reading the paper, A did mental arithmetic on two of the C6 σ conversions (Planck vs full sample, DES vs full sample) and they looked right. C6's expected outcome below is therefore not blind. No code was run and no Ωm was computed before this commit.

## Checks, tolerances and expected outcomes

Paper inputs were DESI v1 Table 1 (rounded means, 2-decimal z_eff). Ours are the DESI DR1 Gaussian likelihood files (unrounded, 3-decimal z_eff). The tolerances allow for that.

| Claim | Check | Tolerance for `held:reproduced` | Expected |
|---|---|---|---|
| C1 | SymPy: re-derive the ratio from Eq. (1); limits Ωm → 0 (ratio = z), Ωm = 1 (ratio = 2[(1+z)^{3/2} − (1+z)]), z → 0 (ratio → z); show d(ratio)/dΩm has one sign on the physical range; find the ratio range reachable for Ωm in [0,1] and in the wider range where E²>0 on [0,z] | Exact | Holds; monotonic decreasing in Ωm |
| C2 | Ratio-method Monte Carlo, 10,000 draws per bin from the 2×2 DR1 covariance; solve Eq. (2) for Ωm by root finding; record the fraction of draws with no root and handle them two ways (drop; clip to range edge) | Central within 0.25 × the paper's error on the facing side; each error bar within 20 % | Reproduces (LRG 0.51 near 0.65) |
| C3 | 50,000 draws, P(Ωm < 0.322); σ one-sided and two-sided | p in [0.011, 0.019]; σ rounds to 2.2 one-sided | Reproduces, one-sided |
| C4 | Per-bin 2D posterior on a dense grid (exact marginals) and an MCMC cross-check; priors as printed. Then widen Ωm prior to [−0.5, 3] (E² > 0 enforced) | As C2 | Reproduces. Widening may move the z = 0.51 upper tail; if central or upper error moves > 0.25σ, the "neither prior impacts our results" sentence is a note (S0) |
| C5 | Full DR1 (12 numbers, full covariance) flat-ΛCDM fit, grid marginal + MLE; with 2-decimal z (rounded from the files) and 3-decimal z; with and without the two z = 0.51 numbers | Central within 0.003, errors within 0.002, MLE within 0.002. Footnote holds if the 3-dp vs 2-dp shift has the right sign and is between 0.0025 and 0.01 | Reproduces |
| C6 | Recompute every σ from printed numbers (convention in conventions.json), then from our C4/C5 numbers | Each σ within ±0.1 of the printed one-decimal value (±0.15 for "~2σ") | Reproduces |
| C7 | Flat-ΛCDM fit to three bins: [BGS 0.295 D_V, LRG 0.51 D_M+D_H], [LRG 0.706, LRG+ELG 0.93], [ELG 1.317, QSO 1.491 D_V, Lyα 2.33] with their covariance blocks; shifts as Δ/sqrt(σ_facing²+σ_facing²) | Centrals within 0.25σ; shifts within ±0.2σ | Reproduces |
| C8 | emcee, flat w0waCDM on DR1 alone, priors as printed, Ωm in [0,1], H0r_d free; with and without the z = 0.51 numbers; 32 walkers, ≥ 20,000 steps after burn-in, R̂/autocorrelation reported | w0 median within 0.10 (with LRG) and 0.15 (without); claim "LRG 0.51 drives w0 > −1" held if removing it moves w0 by ≥ 0.25 towards −1 and −1 lies inside the 68 % interval without it | Reproduces; posterior prior-dominated in wa |
| C9 | χ² of a constant-Ωm fit across all DR1 bins (dof = 12 − 2), PTE and Gaussian σ; compare with the C7 shifts and the 2.8σ DESI-vs-DES number | Held if a within-DESI test also gives ≳ 2.5σ. Otherwise `held:discrepancy` at S0–S1 (wording): the 2.8σ is a difference between two datasets' Ωm, not a test of Ωm varying with redshift (G17, G20) | Narrowing note, S0–S1 |

Severity cap: none of C1–C9 can reach S2 unless a printed number in the abstract (Ωm = 0.67, ~2σ, ~1.8σ) or the headline w0 statement fails its tolerance and the corrected value changes the stated conclusion.

## Fidelity gate (G15) and the DR2 re-test

The DR2 re-test runs only if C2 and C4 pass on DR1. It is a data trace, not verification; its result never changes a claim outcome (G17).

Pre-declared metrics on DESI DR2 (`data/desi-dr2`), same code as DR1:
- **M1.** Ratio-method and grid Ωm at LRG z = 0.510; P(Ωm < 0.322) and its σ. The paper predicts the anomaly "will decrease in statistical significance". `consistent` if σ < 1.5, `tension` 1.5–3, `strong_tension` > 3.
- **M2.** The Table II binning on DR2 (QSO z = 1.484 is now anisotropic: use D_M and D_H). Shift bin 2 → bin 3. The paper predicts an increasing Ωm trend at higher redshift remains.
- **M3.** DR2-alone w0waCDM with and without LRG z = 0.510, same priors: does that bin still move w0 by ≥ 0.25 towards −1?

Literature check (G21) after the runs, not before: the authors' own DR2 follow-up arXiv:2504.04417 (cited at l.282), DESI DR2 (2503.14738) and any comment on this paper.

## Controls

- Known-answer: mock D_M/r_d, D_H/r_d at Ωm = 0.3, H0r_d = 10,000 km/s for each DR1 z; the ratio solver and grid fit must return 0.3 to 1e-6 on noiseless mocks.
- Inversion check: re-insert each solved Ωm into Eq. (2), |residual| < 1e-10 (the paper's own footnote check).
- No author code exists to import (G13). All numbers in the report come from scripts saved in `checks/` with outputs in `checks/out/` (G14).

## Stop rules

Round cap 2 per claim. If C2 or C4 fail and the cause is our input rounding, one addendum may switch to DESI v1 Table 1 rounded inputs; that addendum must be committed before the rerun.

## Addendum A1, 2026-10-07 (after V4, post hoc)

Committed after reading Codex's B report and A-vs-B review, before any A re-run. Everything below is post hoc: it fixes or completes A's first pass, it does not change any tolerance, and the first-pass outputs stay in `checks/out/` unchanged. Re-run outputs go to `checks/out-r2/`.

1. `fit_grid.summarize` quantiles move from an interpolated cumulative sum (biased by about half a grid step) to a trapezoid-integrated CDF. C4, C5, C6, C7 and C9 are re-run. Each grid is checked against one at twice the resolution.
2. C4 gets the per-bin MCMC cross-check and the grid mock recovery this protocol promised and the first pass did not run. Mock recovery is judged on the likelihood maximum (optimiser started from the grid MLE), not on a posterior median.
3. C2 removes the Ωm = 50 inversion cap, reports how many draws would have hit it, implements the promised clip alternative, adds physical-domain [0, 1] drop and clip sensitivities, and checks reinsertion residuals on every draw. C3 adds the exact Gaussian tail P(D_M − R(0.322) D_H < 0) beside the Monte Carlo.
4. C6 applies the ±0.1 criterion to the refits too (the first pass printed them without a verdict), using our own C4 fit for the LRG row.
5. C8 seeds the sampler itself (not only the starting positions), runs 4 independently started ensembles per case, and reports rank-normalised split R-hat and bulk/folded on one fixed walker per ensemble, plus τ. Chains are not committed (size); the summaries are.
6. C9 replaces leave-one-out (not a likelihood-ratio test) with same-data nested tests: (a) one Ωm per Table II group with a common H0r_d, 2 extra parameters; (b) one (Ωm, H0r_d) pair per group, 4 extra parameters. Asymptotic χ² calibration, with the caveat that the groups were chosen after seeing the data.
7. Literature (G21): INSPIRE search of the papers citing 2404.08633 for comments and replies.

Independence note: these re-runs were written after reading B's code and results, so they are not independent of B. They complete A's own evidence; they do not count as a second provider.

## Addendum A2, 2026-10-07 (after Codex's review of A1, post hoc)

Committed before any A2 run. Codex's response review (`model-reports/2026-10-07-codex-B/response-review/RESPONSE-REVIEW.md`) found that A1 left four things open. Tolerances are unchanged.

1. **C8 chains.** Re-run `c8_w0wa.py` with the same seeds. Save every walker thinned by 25 (float32, compressed) plus the full fixed-walker traces used for R-hat in `checks/chains/`, with SHA-256 sums and the Python and package versions. Success criterion, fixed now: max rank-normalised R-hat < 1.01 on every parameter in every case. If a case misses, re-run it once at 40,000 steps. If it still misses, mark it `inconclusive`. Exact reproduction is claimed only for the recorded environment: Codex's Mac run with the same seeds differs by up to 0.008 in w0.
2. **Grid refinement coverage.** A1 promised a 2× check for every grid but ran it only for the five baseline C4 grids. `c10_grid_refine.py` covers the rest: the wide-prior C4 grids, C5 (both z roundings, with and without z = 0.51), C7 groups, and the C9 global grid. Criterion: medians and both errors move by less than 0.1 × the protocol tolerance, i.e. 0.025 × the facing error. The C9 nested optimisations are checked by the spread of χ²min over at least 8 random starts, which must be < 1e-6. C6 reuses C5's grids.
3. **Wording.** C3: the paper's p lies just below our three Monte Carlo estimates, not inside them. C6: the printed summaries and our refit differ in both the central value and the error. The cause is not isolated: 0.2891 rounds to 0.289, so it is not rounding to 0.290. C2: the 1e12 bracket guard is a numerical limit, so "guard exhausted" is not a proof that no root exists.
4. **Literature log.** Record the INSPIRE queries and their hit counts in `literature.md`.
