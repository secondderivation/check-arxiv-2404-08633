"""Seeded multi-ensemble emcee runs and convergence diagnostics (addendum A1, 2026-10-07).
R-hat is the rank-normalised split R-hat of Vehtari et al. (2021), max of bulk and folded, computed on one fixed
walker per independently started ensemble: walkers inside one ensemble are coupled, so they are not independent chains."""
import numpy as np
import emcee
from scipy.stats import norm, rankdata


def _rhat(x):  # x: (chains, draws)
    n = x.shape[1]
    w = x.var(axis=1, ddof=1).mean()
    b = n * x.mean(axis=1).var(ddof=1)
    return np.sqrt(((n - 1) / n * w + b / n) / w)


def _z(x):
    r = rankdata(x.ravel()).reshape(x.shape)
    return norm.ppf((r - 0.375) / (x.size + 0.25))


def rank_split_rhat(chains):
    """chains: (n_chains, n_draws). Returns max(bulk, folded) rank-normalised split R-hat."""
    h = chains.shape[1] // 2
    s = np.concatenate([chains[:, :h], chains[:, -h:]], axis=0)
    return max(_rhat(_z(s)), _rhat(_z(np.abs(s - np.median(s)))))


def run_ensembles(logp, start, ndim, steps, burn, n_ens=4, nwalkers=32, seed=0, vectorize=True, keep=None):
    """start(rng, nwalkers) -> initial positions. Seeds both the starts and emcee's own random state.
    keep: optional dict; if given, it receives per-ensemble post-burn chains for archiving (addendum A2)."""
    pooled, fixed, taus = [], [], []
    for e in range(n_ens):
        rng = np.random.default_rng(seed + e)
        p0 = start(rng, nwalkers)
        s = emcee.EnsembleSampler(nwalkers, ndim, logp, vectorize=vectorize)
        s.random_state = np.random.RandomState(seed + 1000 + e).get_state()
        s.run_mcmc(p0, burn + steps, progress=False)
        ch = s.get_chain(discard=burn)  # (steps, walkers, ndim)
        taus.append(s.get_autocorr_time(discard=burn, quiet=True))
        pooled.append(ch.reshape(-1, ndim)); fixed.append(ch[:, 0, :])
        if keep is not None:
            keep.setdefault('chains', []).append(ch)
    fixed = np.array(fixed)
    rhat = [rank_split_rhat(fixed[:, :, k]) for k in range(ndim)]
    med_spread = np.ptp([np.median(p, axis=0) for p in pooled], axis=0)
    return np.concatenate(pooled), dict(rhat=rhat, tau=np.max(taus, axis=0), steps_over_tau=steps / np.max(taus, axis=0),
                                        ensemble_median_spread=med_spread, seed=seed, ensembles=n_ens, steps=steps, burn=burn)


def fmt_diag(d, names):
    return ' '.join(f'{n}:Rhat={r:.4f},tau={t:.0f},N/tau={k:.0f},ens-median-spread={m:.4f}'
                    for n, r, t, k, m in zip(names, d['rhat'], d['tau'], d['steps_over_tau'], d['ensemble_median_spread']))
