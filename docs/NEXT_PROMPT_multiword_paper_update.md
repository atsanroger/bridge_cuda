# Handoff prompt — update the `multiword-en` LaTeX paper from the QDW/QTW findings

> Paste the block below into a new conversation **in the environment where the
> `multiword-en` LaTeX source lives** (it is NOT on the Bridge WSL box). The
> persistent Claude memory is shared across conversations, so the next agent can
> pull the full detail from `MEMORY.md`; this prompt is the self-contained
> summary in case that environment has different memory.

---

You are helping me update the LaTeX source of my paper (working name **`multiword-en`**)
on an FP32/FP64-base **multiword-precision Domain-Wall Fermion (DWF) solver** in
Bridge++ (CUDA). I have just finished a round of implementation + verification on
the solver side; your job is to revise the paper's results, figures, and framing
to match what is now actually measured and verified. **Do not invent numbers —
ask me for a figure/log if you need the exact value.** First read your persistent
memory index `MEMORY.md` (entries prefixed `qtw_`, `qdw`, `tiled_gauge_caveat`,
`project_prof_meeting_supercomputer`) — it has the primary record. Then propose
edits to the LaTeX before applying them.

## What changed on the solver side (fold these into the paper)

1. **Genuine multiword arithmetic.** The triple-word (TW) path previously capped
   near double accuracy because `tw_add/tw_mul/tw_scal` were only dual-word; they
   are now genuine ε³. The accuracy-critical step is the *device* `tw_mul` of the
   preconditioner coefficient products — precomputing them in host double would
   cap the result at double and must NOT be claimed.

2. **Double-base QTW is finished (new headline result).** A host **triple-double
   (TD, ~159-bit) constant pipeline** feeds the double-build TW preconditioner.
   On a 4⁴×8 true-residual floor: double **FP 1.59e-30, DW 1.03e-62, TW 8.57e-76**
   — three cleanly separated curves, double-triple ≈ **46–60 orders past native
   FP64**, zero NaN. (Figure: `docs/double_base_residual_compare.png`.)

3. **FP32-base results.** TW reaches **5.6e-17, beating native FP64 (3.9e-16),
   FP32-only** (float-float / float-triple). DW gives FP64-level accuracy at
   ~FP64 speed; TW beats FP64 at ~1.66×. Kernels are memory-bound.

4. **The central scientific claim: precision = attainable-accuracy REACH, not
   speed.** Precision sets the **floor**; the spectrum / condition number κ sets
   the **rate**. In our eo-preconditioned DWF the recursive CG residual *overlaps*
   across all six precisions (cross-precision spread <0.18 dex); only the TRUE
   residual ‖b−Ax‖ floor ladders by ε. This was verified on a **REAL thermalized
   16⁴ quenched config** (β=6.0, plaq≈0.596), not only the tiled test gauge — the
   descent overlap is physical (all six precisions hit 1e-12 at the same iter).
   Figures: `docs/multiword_precision_ladder.png`,
   `docs/sweep16_real_q450/ladder_real_mq{0.10,0.05,0.025}.png`,
   `docs/verify1_recursive_vs_true.png`.

5. **Relation to Mukunoki–Ozaki (arXiv:2510.13536).** They run *unpreconditioned*
   CG on ill-conditioned SuiteSparse SPD (κ up to ~1e10) where FP32 loses
   orthogonality, so precision **cuts the iteration count**. We are the opposite
   regime: eo-preconditioning keeps effective κ modest, so precision buys reach,
   not iterations. A toy unpreconditioned-CG κ-sweep isolates the mechanism:
   κ=10 → FP32 and FP64 converge at the same rate; κ≳1e3 → FP32 stalls. Rate-split
   onset is κ ≳ few×10³ (worst case); for real clustered DWF spectra expect
   ~10⁴–10⁵. Figure: `docs/toy_cg_kappa_sweep.png`.

6. **Condition-number anchor from my own α paper (arXiv:2602.14454).** On the
   exact setup (16⁴, β=6.0, eo-preconditioned D†D, double-precision Lanczos):
   at **m=0.001, κ ≈ 1.7e5 (α=0.6 optimal) up to ~1.6e6 (α=0.2), ~3–6e5 at α=1**.
   Our multiword tests sat at m=0.025–0.10 where κ~10³ — two orders below the
   split threshold, which is exactly why all precisions overlapped. Note 16⁴ shows
   finite-volume *saturation* of the mass dependence; a genuinely high-κ, small-mass
   demonstration needs 32⁴+. The α parameter is a κ knob: for the speed/iteration
   story use α≈0.6 (low κ); to **demonstrate the multiword reach/Ozaki advantage**,
   raise κ (α=1 or 0.2 + small mass + large volume).

7. **Honest caveats to keep in the paper.** The standard test gauge is a periodic
   4⁴×8 self-copy (Bloch-decomposable → artificially benign κ), so do NOT claim
   physical volume scaling from it. All quantitative results are on small lattices;
   size scaling is untested and the reconstruction overhead (−43% at this size)
   could change sign at production volume. The surviving *performance* argument is
   double-base QTW on FP64-fast A100/H100 (FP64 1:2), which recovers the
   modest-overhead beyond-FP64 regime; on a 3080 (FP64 1:64) timing is meaningless
   but correctness/precision are demonstrable.

## Task
Walk the `multiword-en` source, then propose (and after I confirm, apply) edits
that (a) add the double-triple/TD result and figure, (b) sharpen the
reach-vs-rate framing with the recursive-vs-true-residual evidence, (c) add the
κ / Ozaki-regime discussion anchored to arXiv:2602.14454, and (d) keep the honest
small-lattice/tiled-gauge caveats. Match the paper's existing notation and section
structure. Show me a diff per section before committing.
