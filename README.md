# Does DESI 2024 Confirm ΛCDM?

[2404.08633v3](https://arxiv.org/abs/2404.08633v3) · Review 0.2.0-rc2

**Reproduced in scope.** The scoped DESI fits reproduce with minor notes and one strict numerical exception.

**Limit:** Some comparisons depend on rounded inputs and priors; the strict C6 tolerance is narrowly missed.

## Checked scope

- Eqs. (1)–(2): the D\_M/D\_H ratio as a function of Ωm, its limits and inversion
- Tables I–III: per-bin Ωm (ratio method and MCMC) and the binned fits
- The p-value and σ for the z = 0.51 bin against Planck
- Full-sample Ωm and the tensions with Planck, Pantheon+, Union3 and DES
- Section III: DESI-alone w0waCDM with and without the z = 0.51 bin, conditional on its printed priors; wa remains weakly constrained
- The 'under pressure at close to 3σ' sentence

Not checked: The survey of other probes (lensing time delays, quasars, GRBs, supernova samples); DESI's own BAO pipeline; The Planck and supernova Ωm values, which are taken as printed

## Review record

| Claim | Method | Outcome / adequacy | Evidence |
|---|---|---|---|
| C1: C1 in fixed reviewed version — Ratio D\_M/D\_H invertible in Ωm | source\_check | held:reproduced / pass | [checks/out/c1.out](checks/out/c1.out) |
| C2: C2 in fixed reviewed version — Table III ratio-method Ωm per bin | source\_check | held:reproduced / pass | [checks/out-r2/c2c3_dr1.out](checks/out-r2/c2c3_dr1.out) |
| C3: C3 in fixed reviewed version — p = 0.0156 (2.154σ one-sided) against printed 0.0147, ≥ 2.2σ; rounding-boundary pass | source\_check | held:reproduced / pass | [checks/out-r2/c2c3_dr1.out](checks/out-r2/c2c3_dr1.out) |
| C4: C4 in fixed reviewed version — N1: Ωm ≤ 1 prior cuts the LRG1 upper tail (+0.18 → +0.23); ~2σ Planck discrepancy survives | source\_check | held:discrepancy / pass | [checks/out-r2/c4.out](checks/out-r2/c4.out) |
| C5: C5 in fixed reviewed version — N2: disclosed rounded z\_eff; release-file z gives DES 2.53σ / 2.50σ, removal effect not robust | source\_check | held:discrepancy / pass | [checks/out-r2/c5_dr1.out](checks/out-r2/c5_dr1.out) |
| C6: C6 in fixed reviewed version — Full DESI vs DES refit 2.802σ against printed 2.7σ, 0.002σ past tolerance; no impact | source\_check | held:discrepancy / pass | [checks/out-r2/c6.out](checks/out-r2/c6.out) |
| C7: C7 in fixed reviewed version — Reproduced | source\_check | held:reproduced / pass | [checks/out-r2/c7_dr1.out](checks/out-r2/c7_dr1.out) |
| C8: C8 in fixed reviewed version — w0wa fit reproduced conditional on printed priors | source\_check | held:reproduced / pass | [checks/out-r3/c8_dr1_2dp.out](checks/out-r3/c8_dr1_2dp.out) |
| C9: C9 in fixed reviewed version — DESI vs DES ~2.8σ on the paper inputs; former N3 kept as context only (cross-dataset, not a DESI-only trend) | source\_check | held:reproduced / pass | [checks/out-r2/c9_dr1.out](checks/out-r2/c9_dr1.out) |
| N4: N4 in fixed reviewed version — l.222 equates D\_M/D\_H with α\_AP; with DESI definitions it is (D\_M,fid/D\_H,fid)/α\_AP. Calculations unaffected | source\_check | held:discrepancy / pass | [EVIDENCE.md](EVIDENCE.md) |

[Protocol and provenance](PROTOCOL.md). This review record is retrospective.

## Review provenance

Prepared with AI assistance. Historical Claude first pass and Codex independent review/critique remain distinct; this assembly reads saved evidence and is neither blind B nor new P2. No scientific code was run in this assembly. Before the replay command below, extract legacy-original.zip into a new directory and read its README.md and requirements.txt for installation, dependency, data and scope instructions.

Tools: Claude (historical A); Codex (historical B/C and current informed assembly). Human involvement: Owner selected preparation scope; no publication decision is recorded for this batch..

## Replay

~~~sh
mkdir replay-work && python3 -m zipfile -e legacy-original.zip replay-work && cd replay-work && python3 reproduce.py --all
~~~

See original source-payload README and recorded P2; not rerun during assembly.. Unpacks the unchanged original candidate into a new work directory and runs its documented scoped replay. Package validation checked the artifact and inherited the earlier replay after verifying the relevant files were unchanged; it did not rerun the scientific calculation.
