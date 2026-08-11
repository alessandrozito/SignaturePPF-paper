#!/usr/bin/env python3
"""
run_tensorsignatures.py -- fit TensorSignatures to the count tensor exported by
R/Comparison_TensorSignatures.R, and write the pieces needed to compare it with
the Poisson Process Factorization (PPF).

ENVIRONMENT
-----------
TensorSignatures 0.5.0 needs tensorflow<=1.15, which has no wheel for Python>=3.8.
Run this with the isolated Python 3.7 interpreter (see the preamble of
R/Comparison_TensorSignatures.R for how it was built):

    $HOME/miniconda3/envs/tensorsig/bin/python run_tensorsignatures.py \
        --tag icgc_H3K36me3 --rank 8

Verified against: tensorsignatures 0.5.0, tensorflow 1.15.0, numpy 1.16.6,
protobuf 3.20.3.

TENSOR LAYOUT (verified against the installed package, not assumed)
------------------------------------------------------------------
TensorSignature.__init__ reads
    samples          = snv.shape[-1]
    p (= 96)         = snv.shape[-2]
    genomic states   = snv.shape[2:-2]
and ts.TensorSignatureData(dimensions=[2]) emits (3, 3, 2, 96, n). So

    snv.shape = (3, 3, n_states, 96, n_samples)
                 ^  ^
                 |  replication strand: 0 leading, 1 lagging, 2 unknown
                 transcription strand:  0 plus,    1 minus,   2 unknown

The R exporter writes 1-based strand indices with all mass in the "unknown"
cell (3, 3) -> here (2, 2), because the mutation channels have already been
pyrimidine-normalised and the strand orientation is not recoverable. That is a
supported mode: index 2 is TensorSignatures' own "unknown" state. The fitted
strand amplitudes (`a`, `b`) are therefore uninformative BY CONSTRUCTION and are
not written out.

WHAT COMES OUT, AND WHY
-----------------------
The comparison quantity is the log intensity ratio between genomic states, per
signature, because that is what both models estimate:

    PPF               beta_kl * (E[x_l | state 2] - E[x_l | state 1])
    TensorSignatures  log( amplitude of state 2 / amplitude of state 1 )

In the fitted object the state amplitude of signature k in state l is
`S[2, 2, l, :, k, 0].sum()`, with state 1 normalised to 1 (equivalently, `k0`
holds the amplitudes of states 2..L relative to state 1; this script derives them
from S so the normalisation is explicit and checkable).

Outputs, all written to the same directory as the input:
    ts_signatures.tsv         96 x rank  spectra (state- and strand-marginal)
    ts_state_amplitudes.tsv   per (signature, state) amplitude and log-ratio
    ts_exposures.tsv          rank x samples
    ts_fit_summary.tsv        rank, log-likelihood, #params, BIC, epochs
    ts_predicted_state_sample.tsv  expected counts per (state, sample); the R side
                              spreads these over the segments to get a per-bin rate
"""

import argparse
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import tensorflow as tf
import tensorsignatures as ts

# Kept in step with config.R: DIR_TENSORSIG. Both honour the same environment
# variable, so moving the project needs no edit to either file.
DEFAULT_DIR = os.path.join(
    os.environ.get("SIGNATUREPPF_PAPER",
                   os.path.expanduser("~/SignaturePPF-paper")),
    "output", "Comparison_TensorSignatures")

# indices of the "unknown" strand state on the two strand axes (0-based)
TX_UNK, REP_UNK = 2, 2


# ----------------------------------------------------------------------
# input
# ----------------------------------------------------------------------
def load_tensor(indir):
    """Rebuild the (3, 3, n_states, 96, n_samples) tensor from the long table."""
    cnt = pd.read_csv(os.path.join(indir, "snv_counts_long.tsv.gz"), sep="\t")
    key = pd.read_csv(os.path.join(indir, "state_key.tsv"), sep="\t")

    samples = sorted(cnt["sample"].unique())
    channels = sorted(pd.unique(cnt["channel"]))
    n_state = int(key["index"].max())

    # Guard against an axis mix-up in the exporter: if `state` and `channel` are
    # swapped upstream, these assertions fire instead of the counts being
    # silently mislabelled.
    if len(channels) != 96:
        raise SystemExit(
            "expected 96 mutation channels, found {} ({}...). The `state` and "
            "`channel` columns are probably swapped in the exported table."
            .format(len(channels), channels[:4]))
    bad = set(cnt["state"].unique()) - set(range(1, n_state + 1))
    if bad:
        raise SystemExit(
            "state column has values {} outside 1..{} declared in state_key.tsv"
            .format(sorted(bad)[:5], n_state))

    s_idx = {s: i for i, s in enumerate(samples)}
    c_idx = {c: i for i, c in enumerate(channels)}

    snv = np.zeros((3, 3, n_state, 96, len(samples)), dtype=np.float64)
    np.add.at(
        snv,
        (cnt["tx_strand"].values - 1,
         cnt["rep_strand"].values - 1,
         cnt["state"].values - 1,
         cnt["channel"].map(c_idx).values,
         cnt["sample"].map(s_idx).values),
        cnt["count"].values,
    )

    total = cnt["count"].sum()
    assert abs(snv.sum() - total) < 1e-6, "counts lost while building the tensor"
    print("snv tensor {}  total counts {:,}".format(snv.shape, int(snv.sum())))
    if snv[TX_UNK, REP_UNK].sum() == snv.sum():
        print("  all mass in the unknown-strand cell: strand axes carry no "
              "information (expected for this pipeline)")
    return snv, samples, channels, key


def load_exposure(indir, snv_shape, samples):
    """Build the normalisation tensor N from state_exposure.tsv.

    N multiplies the predicted counts (`self._Chat1 *= (self.N + 1e-6)`), so it is
    a Poisson/NB OFFSET. It does two jobs here:

      1. copy number. TensorSignatures has no copy-number term, which would make
         the comparison with PPF unfair. state_exposure.tsv holds
         `bin_weight * CN/2` summed per (state, sample) -- the very same exposure
         the PPF model uses -- so supplying it puts both models on equal footing.

      2. masking the uninformative strand cells. `_A[2,2] == 1` is the neutral
         strand state, but the model still predicts counts in all 9 strand cells,
         while our data has mass only in [2,2]. Left alone, the optimiser has to
         drive the strand amplitudes to -inf to explain 8 empty cells. Setting
         N = 0 outside [2,2] makes those predictions ~1e-6 and removes the problem.

    The tensor is scaled to mean 1 over the informative cell: only the RELATIVE
    variation of the offset is identified, the overall level is absorbed by E.
    """
    f = os.path.join(indir, "state_exposure.tsv")
    if not os.path.exists(f):
        raise SystemExit("no state_exposure.tsv in {} -- rerun the R export"
                         .format(indir))
    ex = pd.read_csv(f, sep="\t")
    n_state, n_samp = snv_shape[2], snv_shape[4]
    s_idx = {s: i for i, s in enumerate(samples)}
    M = np.zeros((n_state, n_samp))
    M[ex["state"].values - 1, ex["sample"].map(s_idx).values] = ex["exposure"].values
    if not np.all(M > 0):
        n_zero = int((M <= 0).sum())
        print("  warning: {} (state, sample) cells have zero exposure".format(n_zero))
        M[M <= 0] = M[M > 0].min() * 1e-6
    M = M / M.mean()                                    # relative offset

    N = np.zeros(snv_shape)
    N[TX_UNK, REP_UNK] = M[:, None, :]                  # broadcast over channels
    print("N offset: shape {}  range [{:.3g}, {:.3g}] on the informative cell"
          .format(N.shape, M.min(), M.max()))
    return N


# ----------------------------------------------------------------------
# fitting
# ----------------------------------------------------------------------
def fit_ts(snv, rank, epochs, seed, log_step, N=None, q_other=2, verbose=True):
    """Fit one TensorSignatures model.

    `other` holds the non-SNV mutation types (indels / MNVs / SVs) in the full
    TensorSignatures model. We have none, so a zero matrix is passed: it keeps the
    interface satisfied without contributing to the SNV part of the likelihood.
    q_other must be >= 2 -- the model builds a variable of shape (q - 1, rank).
    """
    other = np.zeros((q_other, snv.shape[-1]))
    model = ts.TensorSignature(
        snv=snv, other=other, rank=rank, N=N,
        objective="nbconst",          # negative binomial, as in the paper
        epochs=epochs, log_step=log_step, display_step=log_step,
        seed=seed, verbose=verbose, id="PPFcmp",
    )
    # Own the session rather than letting fit() create and drop it: the fitted
    # expected counts (Chat1) are a TF tensor, so evaluating them afterwards needs
    # the same session. (fit() documents that it returns the session, but it
    # actually returns the result object -- hence the explicit construction.)
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    model.fit(sess=sess)
    return model, sess


# ----------------------------------------------------------------------
# extraction
# ----------------------------------------------------------------------
def extract(result, samples, channels, key):
    """Pull out the spectra, the genomic-state amplitudes and the exposures."""
    S = np.asarray(result.S)                 # (3, 3, n_state, 96, rank, 1)
    E = np.squeeze(np.asarray(result.E), -1)  # (rank, n_samples)
    n_state, rank = S.shape[2], S.shape[4]
    sig_names = ["TS{:02d}".format(k + 1) for k in range(rank)]

    # --- amplitude of each (signature, state): total mass over the 96 channels
    amp = np.array([[S[TX_UNK, REP_UNK, l, :, k, 0].sum()
                     for k in range(rank)] for l in range(n_state)])   # state x rank
    ref = amp[0]                                                        # state 1 = reference

    rows = []
    for l in range(n_state):
        for k in range(rank):
            rows.append(dict(signature=sig_names[k], state=l + 1,
                             amplitude=amp[l, k],
                             log_ratio_vs_state1=float(np.log(amp[l, k] / ref[k]))))
    # bring in only the readable state name; merging the whole key would collide
    # on its own `state` column and rename ours to state_x/state_y
    if "name" in key.columns:
        lut = key[["index", "name"]].rename(columns={"index": "state"})
    else:
        lut = key.rename(columns={"index": "state"})
        lut = lut[["state"] + [c for c in lut.columns if c != "state"][:1]]
    amp_df = pd.DataFrame(rows).merge(lut, on="state", how="left")

    # --- spectra: marginal over states (and the unknown-strand cell), normalised
    spec = S[TX_UNK, REP_UNK, :, :, :, 0].sum(axis=0)                   # 96 x rank
    spec = spec / spec.sum(axis=0, keepdims=True)
    spec_df = pd.DataFrame(spec, index=channels, columns=sig_names)
    spec_df.index.name = "channel"

    exp_df = pd.DataFrame(E, index=sig_names, columns=samples)
    exp_df.index.name = "signature"
    return spec_df, amp_df, exp_df


def predicted_state_sample(model, sess, samples):
    """Expected counts per (state, sample) under the fitted model.

    Chat1 is TensorSignatures' fitted intensity, shape (3, 3, state, 96, sample).
    Only the [2, 2] strand cell carries signal here, and summing over the 96
    channels gives the per-(state, sample) expected count -- the finest genomic
    resolution TensorSignatures can express, since it has no notion of position
    within a state. The R side spreads these over the ChromHMM segments in
    proportion to each segment's exposure to obtain a per-bin rate comparable to
    PPF's Lambda.
    """
    chat = sess.run(model.Chat1)
    leak = chat.sum() - chat[TX_UNK, REP_UNK].sum()
    if leak > 1e-3 * max(chat.sum(), 1):
        print("  warning: {:.3g} of the predicted mass sits outside the "
              "informative strand cell".format(leak))
    pred = chat[TX_UNK, REP_UNK].sum(axis=1)          # (state, sample)
    rows = [dict(state=l + 1, sample=samples[j], predicted=float(pred[l, j]))
            for l in range(pred.shape[0]) for j in range(pred.shape[1])]
    return pd.DataFrame(rows)


def n_parameters(model, result, rank):
    """Free-parameter count, following tensorsignatures' own definition
    (util.py, Cluster.parameters): 4*95 spectrum parameters plus one parameter
    per strand / genomic-state / mixing / overdispersion dimension, all times the
    rank, plus one exposure per (signature, sample)."""
    d = result.to_dic()
    p = 4 * 95
    for v in ("_a0", "_b0", "_m0", "_T0"):
        p += np.asarray(d[v]).shape[0]
    p += int(np.sum(model.card))          # genomic-state dimensions
    p = p * rank
    p += model.samples * rank
    return int(p)


def fit_summary(model, result, rank):
    """Log-likelihood, parameter count and BIC -- used to choose the rank."""
    logL = float(np.asarray(result.log_L).ravel()[-1]) \
        if np.size(result.log_L) else np.nan
    params = n_parameters(model, result, rank)
    obs = int(model.observations)
    bic = np.log(obs) * params - 2 * logL if np.isfinite(logL) else np.nan
    aic = 2 * params - 2 * logL if np.isfinite(logL) else np.nan
    return pd.DataFrame([dict(rank=rank, log_likelihood=logL, parameters=params,
                              observations=obs, AIC=aic, BIC=bic,
                              epochs=model.epochs, objective=model.objective)])


# ----------------------------------------------------------------------
def run(indir, rank, epochs, seed, log_step, use_exposure=True, outdir=None):
    outdir = outdir or indir
    os.makedirs(outdir, exist_ok=True)
    snv, samples, channels, key = load_tensor(indir)
    N = load_exposure(indir, snv.shape, samples) if use_exposure else None
    if N is None:
        print("  NO exposure offset: the fit is NOT copy-number corrected")
    print("fitting rank {} for {} epochs ...".format(rank, epochs), flush=True)
    model, sess = fit_ts(snv, rank, epochs, seed, log_step, N=N)
    result = model.result

    spec_df, amp_df, exp_df = extract(result, samples, channels, key)
    summ = fit_summary(model, result, rank)

    spec_df.to_csv(os.path.join(outdir, "ts_signatures.tsv"), sep="\t")
    amp_df.to_csv(os.path.join(outdir, "ts_state_amplitudes.tsv"), sep="\t", index=False)
    exp_df.to_csv(os.path.join(outdir, "ts_exposures.tsv"), sep="\t")
    summ.to_csv(os.path.join(outdir, "ts_fit_summary.tsv"), sep="\t", index=False)
    predicted_state_sample(model, sess, samples).to_csv(
        os.path.join(outdir, "ts_predicted_state_sample.tsv"), sep="\t", index=False)
    # keep the run self-contained so the R side (read_ts_fit /
    # compare_chromatin_effects) can point straight at `outdir`
    if os.path.abspath(outdir) != os.path.abspath(indir):
        key.to_csv(os.path.join(outdir, "state_key.tsv"), sep="\t", index=False)

    print("\n--- state amplitudes (log ratio vs state 1) ---")
    print(amp_df[amp_df.state > 1][["signature", "state", "log_ratio_vs_state1"]]
          .round(3).to_string(index=False))
    print("\n--- fit ---")
    print(summ.round(1).to_string(index=False))
    print("\nwritten to {}".format(outdir))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True,
                    help="sub-directory written by export_ts_input(), "
                         "e.g. icgc_H3K36me3")
    ap.add_argument("--dir", default=DEFAULT_DIR, help="parent output directory")
    ap.add_argument("--rank", type=int, default=8, help="number of signatures")
    ap.add_argument("--epochs", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None,
                    help="directory for the outputs (default: the input "
                         "directory). Use a per-rank sub-directory when sweeping "
                         "the rank, otherwise each fit overwrites the previous one")
    ap.add_argument("--no-exposure", action="store_true",
                    help="fit WITHOUT the copy-number / strand-mask offset N "
                         "(not recommended: the comparison with PPF is then "
                         "not copy-number corrected)")
    ap.add_argument("--log-step", type=int, default=100,
                    help="must be <= epochs, otherwise the log stays empty and "
                         "the fit errors out when it is read")
    args = ap.parse_args()

    indir = os.path.join(os.path.expanduser(args.dir), args.tag)
    if not os.path.isdir(indir):
        raise SystemExit("no such directory: {}\nRun the R export step first."
                         .format(indir))
    if args.log_step > args.epochs:
        raise SystemExit("--log-step must be <= --epochs")
    run(indir, args.rank, args.epochs, args.seed, args.log_step,
        use_exposure=not args.no_exposure,
        outdir=os.path.expanduser(args.out) if args.out else None)


if __name__ == "__main__":
    main()
