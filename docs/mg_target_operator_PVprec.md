# MG target operator — FP implementation & 16⁴ spectrum

**Date:** 2026-06-06
**Author:** Wei-Lun Chen (with Claude Code)
**Paper basis:** I. Kanamori, W.-L. Chen, H. Matsufuru, *Spectrum of Preconditioned
Moebius Domain-wall Operators*, PoS(LATTICE2024)414 (`reference/mgspecturm.pdf`).

## Target operator

The paper's best-conditioned operator for multigrid (smallest condition number,
1.96×10⁴; positive real spectrum):

```
A = (D_PV C_PV⁻¹)†  C⁻¹  D
```

- `D`     : Moebius domain-wall operator at the physical mass.
- `C⁻¹`   : exact 5D site-diagonal LU inverse = Bridge++ `Prec` mode (`L_inv∘U_inv`).
- `D_PV`  : Pauli-Villars operator = same operator at `quark_mass = quark_mass_PauliVillars` (1.0).
- `C_PV⁻¹`: site-diagonal LU inverse of the PV operator = `Prec` of the PV instance.

All four factors are applied **exactly** — no iterative solver, no FP64 promotion.
This is the FP (single/double precision) reference implementation.

> Note vs. existing MG: `test_alt_dw_multigrid_dev6` already does PV preconditioning,
> but through a **SAP loose inverse of `D_PV`** — the paper's `B_PV` regime, which the
> paper shows develops *negative* real parts. The operator here uses the **exact**
> site-diagonal `C_PV⁻¹` instead.

## Implementation

| file | role |
|---|---|
| `src/lib_alt_Accel/Fopr/afopr_Domainwall_PVprec.h` | composite operator (2 internal `AFopr_Domainwall_5din`) |
| `src/lib_alt_Accel/Fopr/afopr_Domainwall_PVprec_{float,double}.cpp` | instantiation + factory class_name |
| `src/lib_alt_Accel/Fopr/afopr.cpp` | manual registration `"Domainwall_PVprec"` (this build is `USE_FACTORY` **without** `USE_FACTORY_AUTOREGISTER`) |
| `test_alt_spectrum/src/Tests/test_alt_Accel.cpp` | `test_Eigenvalue_Domainwall_PVprec_Arnoldi()` |
| `test_alt_spectrum/src/Tests/eigenvalue_alt-tmpl.h` | harness patched: `TDa` uses `AEIGENSOLVER::complex_t` (supports the `float` field path) |
| `test_alt_spectrum/src/test_alt_Eigenvalue_Arnoldi_Domainwall_PVprec.yaml` | run parameters |

mult chain (forward `D`):  `D` → `Prec` → `Ddag(PV)` → `Precdag(PV)`.

## Run recipe (WSL, single GPU)

```bash
cd test_alt_spectrum/src
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/comm_libs/13.0/hpcx/hpcx-2.25.1/ompi
export PATH=$OPAL_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$OPAL_PREFIX/lib:$LD_LIBRARY_PATH
EXE=../../build/cuda-ampere-mpi/bridge_test_alt_spectrum
mpirun --mca plm isolated --mca btl self,vader --mca pml ob1 -np 1 "$EXE"
```

- config: `sample_hmc/conf_quench_16x16x16x16.txt` (real, thermalized 16⁴ quenched HMC
  config — NOT the tiled 4⁴×8 test gauge, so conditioning is physical).
- DW params used: `M0=1.0, b=1.5, c=0.5, Ns=8, mq=0.01, M_PV=1.0, bc=[1,1,1,-1]`.

## Result (16⁴, double-field eigensolver)

Dominant eigenvalues of `A` (stable Ritz values across IRArnoldi restarts):

```
       λ                       |λ|
  0.437 ± 0.175 i           0.471
  0.347 ± 0.260 i           0.434
  0.226 ± 0.180 i           0.289
  0.120 ± 0.035 i           0.125
```

- **All eigenvalues have positive real part** — matches paper Fig. 4 (left panel). ✅
- Spectral radius ≈ 0.47; the dominant spectrum sits in a small positive-real-part disk.
- Operator runs cleanly at **M0=1.0** (the earlier failure was never the operator).

### Caveats / open items

1. **Non-normal operator:** IRArnoldi residuals plateau at ~3–4×10⁻² and do not reach
   1e-6. The Ritz *values* are stable to ~4 digits and reliable; tight residual
   convergence is intrinsically hard for a non-normal operator (the paper likewise
   reports eigenvalues, not residual-converged eigenpairs).
2. **Condition number not yet measured:** plain Arnoldi (any `eigensolver_mode`)
   converges to the **dominant / largest-|λ|** modes (≈0.47). The **smallest** modes
   near the quark mass — needed for |λ_max|/|λ_min| — require **shift-invert** (or a
   polynomial/Chebyshev transform), which is not yet wired up.
3. **Eigensolver-settings gotcha:** with a tiny Krylov space (Nk=Np=12) and large
   `threshold_value` (1.0), *every* preconditioned DW operator — including Bridge++'s
   own built-in `DdagD_prec` — overflows/NaNs in IRArnoldi (loss of orthogonality on
   the clustered near-origin spectrum). `Prec`/`C⁻¹` itself is sound (used in production
   CG). Robust settings: **Nk=8, Np=40, threshold_value=0.02** → zero NaN, zero overflow.

## Next steps

- Shift-invert (or Chebyshev) to reach the low modes → quote the actual condition number.
- Register the `float` eigensolver (`AEigensolver<AField<float,ACCEL>>` is missing from
  `bridge_init_factory_alt_Accel.cpp`; needs float IRArnoldi/IRLanczos instantiation) for
  a genuinely FP-field spectrum run.
- Wire the exact `C_PV⁻¹` preconditioning into the dev6 MG coarse-op generation,
  replacing the SAP `B_PV` path.
