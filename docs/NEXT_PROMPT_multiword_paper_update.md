# Handoff prompt — update `docs/latex/multiword.tex` from the measured QDW/QTW results

> Paste the block below into a fresh conversation **in this same Bridge repo**
> (the LaTeX source and `pdflatex` are both on this WSL box). The persistent
> The persistent agent memory is shared, so the next agent has the full record; this prompt is
> the targeted brief.

---

I want you to update my LaTeX paper at **`docs/latex/multiword.tex`** (single
file, ~1810 lines, `article` class, builds with `pdflatex docs/latex/multiword.tex`).
The paper is currently a **review + methods + proposal**: it explains multi-word
(DW/TW/QDW/QTW) arithmetic, Wilson/DWF operators, and eo-preconditioning, then in
`\subsection{Future Research...}` (\label{sec:future}) proposes a 2×3 benchmark
with *predicted* hypotheses H1/H2/H3. I have since **implemented and measured** the
DWF QDW/QTW solver, so parts of those predictions now have empirical answers — some
confirming, **some correcting** the paper. Update it carefully and honestly.

**Process (do this first, in order):**
1. Read your persistent memory `MEMORY.md` and the entries prefixed `qtw_`,
   `qdw`, `tiled_gauge_caveat`, `project_prof_meeting_supercomputer` — the primary
   record with exact numbers and caveats.
2. Read `docs/latex/multiword.tex` end to end; note its notation
   ($\uu$, $\kappa(D^\dagger D)$, QDW/QTW, the §Condition Number Analysis threshold
   table, and the H1/H2/H3 hypotheses).
3. Propose a **section-by-section diff and a plan** before editing. Wait for my OK.
4. **Do not invent numbers.** Exact floors/logs are in `docs/sweep16_real_q450/`
   and the memory; ask me or read those if unsure.
5. After editing, run `pdflatex` and confirm it builds clean.

**Reuse existing citation keys:** `Mukunoki2025` (2510.13536), `Chen2026` (my
α/modified-DWF paper, 2602.14454), `Kanamori2026` (FP16), `Chen2025`, `Akahoshi2022`.

## The measured findings, mapped to where they belong in the paper

1. **New empirical Results (add a subsection in §"Multi-Word Arithmetic in Lattice
   QCD Solvers", before §Future Research).** Genuine multiword arithmetic now
   verified end-to-end in Bridge++ CUDA:
   - FP32-base: **TW reaches 5.6e-17, beating native FP64 (3.9e-16), FP32-only**;
     DW ≈ FP64 accuracy at ≈FP64 speed; TW beats FP64 at ~1.66×.
   - **Double-base (new): triple-double (TD, ~159-bit) constant pipeline → double-TW
     floor 8.57e-76**, i.e. three cleanly separated FP/DW/TW curves (double FP
     1.59e-30, DW 1.03e-62, TW 8.57e-76), ~46–60 orders past native FP64, zero NaN.
   - Figures to `\includegraphics`: `docs/double_base_residual_compare.png`,
     `docs/multiword_precision_ladder.png`.

2. **Correct the framing of §"Effective Precision Thresholds" / §Future Research.**
   The paper's threshold theory ($\uu_{\mathrm{eff}}\,\kappa \gtrsim 1$ → stagnation)
   is right *in principle*, but its predicted **iteration-count reduction (H1/H2)
   was NOT observed in our tests** — because at our tested masses the
   eo-preconditioned $\kappa(D^\dagger D)$ was far below the stagnation threshold.
   What we measured instead is **attainable-accuracy REACH**: precision sets the
   residual *floor*; the spectrum/$\kappa$ sets the *rate*. The recursive CG
   residual overlaps across all six precisions (spread <0.18 dex); only the TRUE
   residual $\|b-Ax\|$ floor ladders by $\uu$. Verified on a **real thermalized
   16⁴ quenched config** (not just the tiled test gauge). Reframe H1/H2 from
   "predicted speedup" to "regime-dependent: confirmed reach; iteration savings
   require large $\kappa$ (small mass + large volume), still open."
   Figures: `docs/verify1_recursive_vs_true.png`,
   `docs/sweep16_real_q450/ladder_real_mq{0.10,0.05,0.025}.png`.

3. **Ground §"Condition Number Analysis" in real numbers via `\cite{Chen2026}`.**
   The current κ table is generic ($\kappa(D)=10^5\!-\!10^{13}$). Add the *measured*
   eo-preconditioned DWF values from Chen2026 on the exact setup (16⁴, β=6.0,
   double-Lanczos): at **m=0.001, $\kappa(D^\dagger D)\approx 1.7\times10^5$ (α=0.6
   optimal) up to $\sim1.6\times10^6$ (α=0.2), $\sim3\!-\!6\times10^5$ at α=1**. Note
   our multiword tests sat at m=0.025–0.10 where $\kappa\sim10^3$ — two orders below
   the split threshold, explaining the overlap. Note 16⁴ shows finite-volume
   *saturation* of the mass dependence (need 32⁴+ to push $\kappa$), and that the α
   parameter is a $\kappa$ knob (α≈0.6 minimizes $\kappa$ for speed; α=1 or 0.2
   raises it to expose the multiword/Ozaki regime).

4. **Sharpen the Mukunoki–Ozaki contrast (§"Why Lattice QCD Differs").** They:
   unpreconditioned CG on ill-conditioned SuiteSparse SPD ($\kappa\!\sim\!10^{10}$)
   where FP32 loses orthogonality → precision cuts iterations. Us: eo-preconditioned
   DWF, modest effective $\kappa$ → precision buys reach, not iterations. A toy
   unpreconditioned-CG $\kappa$-sweep isolates the mechanism ($\kappa=10$: FP32=FP64
   rate; $\kappa\gtrsim10^3$: FP32 stalls; onset $\kappa\gtrsim$ few$\times10^3$,
   real clustered DWF likely $10^4\!-\!10^5$). Figure: `docs/toy_cg_kappa_sweep.png`.

5. **Keep / add honest caveats.** The standard test gauge is a periodic 4⁴×8
   self-copy (Bloch-decomposable, artificially benign $\kappa$) → do not claim
   physical volume scaling from it. All results are on small lattices; size scaling
   is untested and the reconstruction overhead (−43% here) may change sign at
   production volume. The surviving *performance* argument is double-base QTW on
   FP64-fast A100/H100 (FP64 1:2); on the 3080 (FP64 1:64) timing is meaningless but
   correctness/precision are demonstrable.

6. **Collapse `\subsection{Memory Layout for GPU Implementation}`
   (\label{sec:layout}, ~lines 1360–1547).** It is currently written in
   design-exploration voice — enumerating "Three Options for DW Extension"
   (A: ri 2→4, B: separate hi/lo SoA, C: double4) with pros/cons and a "Layout
   Comparison" table — but the implementation is now **settled**, so the
   alternatives-to-evaluate framing contradicts the measured-results narrative.
   - **Delete** the Option A/B/C enumeration and the Layout Comparison table.
   - Replace with one short paragraph stating the *implemented* layout: **double4
     `{r_hi,i_hi,r_lo,i_lo}` for QDW, and SoA×3 (three `double2` hi/mid/lo buffers)
     for QTW**; both keep all words coalesced (32×32 B = 1024 B / 32×16 B = 512 B
     per warp); the 48-byte packed struct is rejected for alignment. Keep the
     `index_d4` equation (\eqref{eq:d4_index}).
   - The two CUDA listings (`qdw_hopping_kernel`, `normalize_*`) are **simplified
     schematics** (e.g. `gauge[/*neighbor index*/0]` is a placeholder), not the
     real kernels — either delete them, move to an appendix clearly labelled
     "schematic", or replace with a real excerpt from
     `src/lib_alt_Accel/BridgeCUDA/`. Do not present toy kernels as the actual
     implementation. Net target: ~190 lines → ~15–20 lines.

## Deliverable
A revised `multiword.tex` that turns the relevant proposal/prediction text into a
*measured-results* narrative — adding the double-triple reach result and figures,
reframing the reach-vs-rate distinction, grounding $\kappa$ in `\cite{Chen2026}`,
sharpening the Ozaki-regime contrast, and preserving the honest small-lattice
caveats — matching the paper's existing notation and section structure, building
clean under `pdflatex`. Show me a per-section diff before applying.
