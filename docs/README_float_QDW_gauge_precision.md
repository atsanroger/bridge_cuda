# Float-base QDW: vector precision is not enough — the gauge field caps accuracy

**Status:** finding documented 2026-05-25. "Before" baseline (float gauge) recorded
here for a future before/after comparison once a float-float gauge is implemented.

## TL;DR (中文)

float-base QDW(float-float,~48-bit)只把**向量**做成 double-double 風格,
但 **gauge link 仍是單精度 float**。hopping 項每次 `qdw_mult_uc(float2 gauge, float4 vec)`
都把 gauge 捨入到 float,所以**算符 A 本身只有 float 精度**,向量再精準也救不回來。
結果:`float_eo + DW` 與 `float_eo + FP` 收斂與結果幾乎完全相同(誤差 ~1e-7),
QDW 向量精度毫無幫助。要讓 float-base 真正模擬 double,**gauge(以及 source/參考)也必須
做成 float-float**。這是 roadmap #3 真正剩下的工作。

## Setup

- Code: Bridge++ alt_Accel (BridgeCUDA), 5D Domain-Wall fermion, even-odd preconditioned.
- GPU: RTX 3080 (consumer; FP64 = 1/64 rate — timings are not meaningful, only correctness).
- Lattice: 32×8×8×12, single rank.
- Solver: CGNR, `convergence_criterion_squared = 1e-28`, `max_iter = 100`, `max_restart = 40`.
- Test: `Spectrum_Domainwall_alt.hadron_2ptFunction`
  (`test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml`).
- Reference ("expected") = double-precision corelib result = **4.58304049059966e-01**.
- `multiword_mode`: `FP` (plain), `DW` (double-word / QDW), `TW` (triple-word reductions only).
- Field `nin`: 192 = single-word field; 384 = double-word (QDW) field.

The solver's "diff" is the **recursive** relative residual² (`|r|²/|b|²` from the CG
recurrence); the test's "diff2" is the **true** residual² recomputed by the
double-precision corelib operator: `diff2 = ‖D·x − b‖² / ‖b‖²`.

## Results ("before" baseline — float gauge everywhere)

| field_type | mw | field nin | result | \|err\| vs double | true diff2 | recursive resid² @conv | test |
|---|---|---|---|---|---|---|---|
| double_eo | FP | 192 | 4.58304049059965e-01 | ~1e-15 | 6.8e-28 | 3.8e-30 | ✓ |
| double_eo | TW | 192 | 4.58304049059965e-01 | ~1e-15 | 6.8e-28 | 3.8e-30 | ✓ |
| double_eo | DW | 384 | 4.58304049059965e-01 | ~1e-15 | 6.4e-28 | 3.7e-30 | ✓ |
| float_eo  | FP | 192 | 4.58304088895693e-01 | 8.7e-08 | 2.8e-13 | 9.410e-15 | ✗ |
| float_eo  | DW | 384 | 4.58304090010420e-01 | 8.9e-08 | 2.3e-13 | 9.401e-15 | ✗ |

(`float_eo + DW` previously crashed — the float QDW kernels were `exit(1)` stubs.
They are now implemented as float-float, see "Implementation status" below.)

### The decisive observation

The per-iteration recursive-residual traces of `float_eo+FP` and `float_eo+DW` are
**near-identical** — e.g. at iteration 202 of the first solve:

```
float_eo + FP :  202   9.410241276529632e-15
float_eo + DW :  202   9.401143295643174e-15
```

If float-float **vector** arithmetic were buying precision, DW would pull away from
FP. It does not. (Note the displayed value is `|r|²/|b|²`, so 9.4e-15 ≈ a *float-level*
relative residual of ~3e-7.)

## Root cause

The QDW even-odd hopping kernels read each gauge link as a single-precision
`real2` (= `float2` in the float build) and multiply it into the float-float vector:

```c
// src/lib_alt_Accel/BridgeCUDA/src/mult_Domainwall_5din_eo_cuda_qdw-inc.h  (DWF_GMUL_FWD/BCK)
real2 _u; real4 _t;
_u.x = u_ptr[IDX2_G_R(0,0,isg)];   // gauge link rounded to float
_u.y = u_ptr[IDX2_G_I(0,0,isg)];
... qdw_mult_uc(_u, vt1_c0);        // float2 (single) × float4 (float-float)
```

So every hopping multiply discards the link's precision below ~1e-7. The operator
`A_float` is therefore only a *float-accurate approximation* of the true Dirac
operator. CG faithfully solves `A_float x = b` to ~float relative residual, but that
`x` is ~1e-7 away from the double answer — and the float-float **vectors** cannot
recover precision the **gauge multiply** already destroyed.

This generalizes: **the achievable accuracy of an iterative solve is bounded by the
precision of the operator's input data (gauge links, source), not by the precision of
the working vectors alone.** Extended-precision vectors only remove *accumulation*
error (round-off stagnation), which is irrelevant when the operator itself is the
floor.

## The improvement to test in the future ("after")

To make float-base actually emulate double, extend precision on the **operator inputs**:

1. **Float-float gauge** (primary): store links as `float4` (hi/lo) and replace the
   `qdw_mult_uc(float2, float4)` link multiply in `DWF_GMUL_FWD` / `DWF_GMUL_BCK`
   with a float-float × float-float complex multiply. Requires a float-float gauge
   field + gauge convert/split that keeps the low word.
2. (Secondary) source and the final reconstruction in float-float, if the observable
   needs it.

### How to compare before vs after

Re-run the exact `float_eo + DW` config below and compare against this table:

- **Before (now):** result ≈ 4.5830409e-01, error ~9e-8, true `diff2` ~2e-13.
- **After (float-float gauge), expected:** result error should drop toward the
  float-float floor (~1e-12…1e-13) and `diff2` toward ~1e-26, approaching the
  `double_eo` rows above. The `float_eo+DW` vs `float_eo+FP` residual traces should
  then **diverge** (DW pulling lower), unlike now.

## Reproduction

In `build/<cfg>/` (e.g. `build/cuda-ampere-mpi/`), edit
`test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml`:

```yaml
Test_Spectrum:
  field_type: float_eo        # or double_eo
Solver:
  solver_type: CGNR
  multiword_mode: DW          # FP | DW | TW
  convergence_criterion_squared: 1.0e-28
  maximum_number_of_iteration: 100
  maximum_number_of_restart: 40
  verbose_level: Detailed     # to dump the per-iteration recursive residual
```

Run (note: use the hpcx mpirun shipped with nvhpc, the plain `/bin/mpirun` fails
`MPI_Init` with the hpcx-linked binary):

```
/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi/bin/mpirun \
  -np 1 ./bridge_test_alt_spectrum
```

The per-iteration `<iter>  <|r|²/|b|²>` lines give the recursive-residual trace; the
final `result/expected/diff` lines give the observable accuracy.

## Update 2026-05-25 — extended_gauge switch built; the bottleneck moves, not vanishes

A runtime switch `extended_gauge` (YAML, `Fopr` block; default false) was added. When
true it carries a **low word** through the gauge links, the source convert, and the
result reverse (so a float base keeps ~double precision there). Implementation:
`qdw_mult_cc` (DD×DD complex) in `bridgeACC_AField.h`; `DWF_GMUL_*_FF` + a
`template<bool EXT>` hopb kernel; `convert_gauge_dd`; `m_Ueo_lo` + flag plumbing in
the eo fopr; lo-word handling in the QDW `convert`/`reverse`.

**Result: `float_eo + DW + extended_gauge=true` is STILL ~1e-7** (result
4.5830409…, error ~1e-7) — essentially unchanged. Reason: the **5D diagonal
coefficients are still single precision.** In the 5dir / LU kernels the Möbius/DWF
coefficients are computed in `real_t` and multiplied into the DD vector as a
single-word factor:

```c
const real_t FF1 = b_con[is]*(4.0 - M0) + 1.0;   // float for a float base
const real_t Fup = ... 0.5*FF2*alpha;            // float
... dw_scal(f1, vt.x, vt.z, ...)                 // single-word coeff × DD vector
```

So the operator's 5D-diagonal entries carry ~1e-7 error regardless of gauge/source
precision, and that now dominates.

**Refined finding (the paper point):** making a float base emulate double is not a
single fix — it requires extended precision on *every* single-precision input to the
operator, removed one layer at a time:

1. Vectors (float-float field) — no effect alone (operator is the floor).
2. + Gauge + source + reverse (the `extended_gauge` switch) — still ~1e-7,
   because the 5D coefficients remain single precision.
3. + 5D coefficients (b,c,mq,M0,alpha and derived dp,dm,dpinv,e,f in float-float,
   with DD×DD coefficient multiplies in the 5dir/LU kernels) — **not yet done**;
   this is the remaining cap.

Each layer only exposes the next floor. The `extended_gauge` switch is the
infrastructure for layers (1)+(2); finishing the emulation needs layer (3).

## Implementation status (float-base QDW vectors)

Implemented 2026-05-25 (compiles + runs; vectors are float-float, gauge still float):
- Templated the double-double primitives (`TwoSum/TwoProd/dw_add/dw_scal/qdw_mult_uc`)
  in `bridgeACC_AField.h` (overloaded `fma`; `dw_scal` takes a separate multiplier
  type so a `double` literal works on a `float` field).
- Parameterized the EO QDW operator kernels
  (`mult_Domainwall_5din_eo_cuda_qdw-inc.h`, `..._eo_inv_cuda_qdw-inc.h`) to
  `real4`/`real2`/`real_t`; added `typedef double4/float4 real4` (+`real2`) per `.cu`;
  `bridgeACC_Domainwall_5din_float.cu` now includes those kernels instead of the
  `exit(1)` stubs. The double path is unchanged (semantically identical for `double`).
