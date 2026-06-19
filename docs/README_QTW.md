# Quasi Triple-Word (QTW) Arithmetic Enhancements for Bridge++

This branch implements **Quasi Triple-Word (QTW)** arithmetic to improve the numerical precision of reduction operations (dot products, norms) in CUDA kernels. This is particularly effective for achieving deep convergence in solvers and eliminating numerical noise in complex correlator measurements.

## Key Improvements

### 1. Robust QTW Core Functions
Added high-performance Error-Free Transformations (EFT) in `afield_cuda-inc.h`:
- `TwoSum(a, b, x, y)`: Computes $x = a + b$ and the exact rounding error $y$.
- `TwoProdFMA(a, b, x, y)`: Computes $x = a \times b$ and the exact rounding error $y$ using FMA.

### 2. Enhanced Reduction Kernels
Implemented QTW versions for all core BLAS reduction operations:
- **`norm2_reduce_fused_kernel_qtw`**: High-precision vector norm.
- **`dot_reduce_fused_kernel_qtw`**: High-precision real dot product.
- **`dotc_reduce_fused_kernel_qtw`**: High-precision complex dot product (critical for physics measurements).

### 3. Precision-Preserving Reduction Path
- **Shared Memory EFT**: The intra-block reduction phase now uses `TwoSum` throughout, preserving error terms until the final aggregation.
- **Dual-Word Output**: Kernels now output high and low precision parts separately for multi-block reductions.
- **Second-Pass QTW**: Added `reduce_kernel_multiblocks_qtw` to aggregate partial results using compensated summation, preventing precision loss in large-scale reductions.

### 4. Runtime Configuration
Introduced a flexible runtime switch:
- **YAML control**: Add `use_qtw: true` under the `Solver` block in your configuration file.
- **Layered Propagation**: The flag is seamlessly propagated from `Parameters` -> `ASolver` -> `AField` -> `BridgeACC` -> CUDA Kernels.

## Verification & Impact

### Numerical Precision
Verification on Domain-Wall solver tests showed:
- **Imaginary Noise Suppression**: Complex correlators (which should be real) saw imaginary components drop from $10^{-20}$ (standard Double) to ~ $10^{-36}$ (effective QTW precision).
- **Solver Stability**: Enabled stable convergence down to $10^{-32}$ and beyond, which is typically impossible with standard double precision due to rounding stagnation.

### Performance Note
On consumer-grade GPUs (like RTX 3080), QTW operations are **Compute-bound** due to the limited FP64 throughput. While QTW ensures maximum accuracy, it may increase execution time per iteration on these cards. Performance gains via reduced iteration counts are most prominent on data-center GPUs (A100/H100) or in ill-conditioned systems.

## Usage
To enable QTW in your solver:
```yaml
Solver:
  solver_type: CGNR
  use_qtw: true
  convergence_criterion_squared: 1.0e-32
```

##
執行方法，在MAIN內make -j 24之後去test_alt_spectrum/src 後再make -j 24
之後到test_alt_spectrum/tests內執行
mpirun -np 1 ./bridge.elf
如果通過測驗的話代表計算正確