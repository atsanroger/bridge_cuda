/*!
      @file    afield_cuda-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#include <stdio.h>

//====================================================================
void afield_init(real_t *data, const int size)
{
  enter_data_create(data, size);
  set(data, 0.0, size);
}

//====================================================================
void afield_tidyup(real_t *data, const int size)
{
  exit_data_delete(data, size);
}

//====================================================================
void copy_to_device(real_t *v, int nv)
{
  update_device(v, 0, nv);
}

//====================================================================
void copy_to_device(real_t *v, int nv1, int nv)
{
  update_device(v, nv1, nv);
}

//====================================================================
void copy_from_device(real_t *v, int nv)
{
  update_host(v, 0, nv);
}

//====================================================================
void copy_from_device(real_t *v, int nv1, int nv)
{
  update_host(v, nv1, nv);
}

//====================================================================
__global__ void set_kernel(real_t *v, real_t a)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  v[i] = a;
}

//====================================================================
void set(real_t *y, real_t a, int size)
{
  real_t *y_dev = (real_t *)dev_ptr(y);

  int nth = VECTOR_LENGTH;
  int nbl = size / nth;

  set_kernel<<<nbl, nth>>>(y_dev, a);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void afield_set_kernel(real_t *v, real_t a, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in = 0; in < nin; ++in)
  {
    v[IDX2(nin, in, site)] = a;
  }
}

//====================================================================
void afield_set(real_t *y, real_t a, const int nin, const int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  afield_set_kernel<<<nbl, nth>>>(y_dev, a, nin);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void convert_kernel(real_t *v, double *w, int nin, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  if (site < nvol)
  {
    for (int in = 0; in < nin; ++in)
    {
      v[IDX2(nin, in, site)] = w[IDX_CORE(nin, in, site)];
    }
  }
  else
  {
    for (int in = 0; in < nin; ++in)
    {
      v[IDX2(nin, in, site)] = 0.0;
    }
  }
}

//====================================================================
void convert(real_t *y, double *x, int nin, int nvol, int nvol_pad)
{
  real_t *y_dev = (real_t *)dev_ptr(y);

  int size2 = nin * nvol;
  enter_data_create(x, size2);

  double *x_dev = (double *)dev_ptr(x);
  update_device(x, 0, size2);

  int nth = VECTOR_LENGTH;
  int nbl = nvol_pad / nth;

  convert_kernel<<<nbl, nth>>>(y_dev, x_dev, nin, nvol);

  cudaDeviceSynchronize();

  exit_data_delete(x, size2);
}

//====================================================================
__global__ void reverse_kernel(double *v, real_t *w, int nin, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  if (site < nvol)
  {
    for (int in = 0; in < nin; ++in)
    {
      v[IDX_CORE(nin, in, site)] = w[IDX2(nin, in, site)];
    }
  }
}

//====================================================================
void reverse(double *y, real_t *x, int nin, int nvol, int nvol_pad)
{
  int size2 = nin * nvol;
  enter_data_create(y, size2);
  double *y_dev = (double *)dev_ptr(y);

  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol_pad / nth;

  reverse_kernel<<<nbl, nth>>>(y_dev, x_dev, nin, nvol);

  cudaDeviceSynchronize();

  update_host(y, 0, size2);

  exit_data_delete(y, size2);
}

//====================================================================
// General FP→QDW conversion for fields with nin4 complex components per site.
// Works for both 4D (nin4=NC*ND=12) and 5D EO (nin4=NC*ND*Ns=96) fields.
// FP layout:  IDX2(2*nin4, 2*in4+{0,1}, site)  for {re,im}
// QDW layout: 4*IDX2(nin4, in4, site)+{0,1,2,3} for {hi_re,hi_im,lo_re,lo_im}
__global__ void copy_to_qdw_kernel(real_t *qdw, const real_t *std_v, int nin4, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site >= nvol) return;

  for (int in4 = 0; in4 < nin4; ++in4) {
    int base  = 4 * IDX2(nin4,    in4,     site);
    int r_idx =     IDX2(2*nin4,  2*in4,   site);
    int i_idx =     IDX2(2*nin4,  2*in4+1, site);
    qdw[base + 0] = std_v[r_idx];
    qdw[base + 1] = std_v[i_idx];
    qdw[base + 2] = real_t(0);
    qdw[base + 3] = real_t(0);
  }
}

//====================================================================
// General QDW→FP conversion (inverse of copy_to_qdw_kernel).
__global__ void copy_from_qdw_kernel(real_t *std_v, const real_t *qdw, int nin4, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site >= nvol) return;

  for (int in4 = 0; in4 < nin4; ++in4) {
    int base  = 4 * IDX2(nin4,   in4,     site);
    int r_idx =     IDX2(2*nin4, 2*in4,   site);
    int i_idx =     IDX2(2*nin4, 2*in4+1, site);
    std_v[r_idx] = qdw[base + 0] + qdw[base + 2];  // hi_re + lo_re
    std_v[i_idx] = qdw[base + 1] + qdw[base + 3];  // hi_im + lo_im
  }
}

//====================================================================
void copy_to_qdw(real_t* qdw_v, const real_t* std_v, int nvol, int nin4)
{
  real_t* qdw_dev = (real_t*)dev_ptr(qdw_v);
  real_t* std_dev = (real_t*)dev_ptr(const_cast<real_t*>(std_v));

  int nth = VECTOR_LENGTH;
  int nbl = CEIL_NWP(nvol) / nth;
  if (nbl == 0) nbl = 1;

  copy_to_qdw_kernel<<<nbl, nth>>>(qdw_dev, std_dev, nin4, nvol);
  cudaDeviceSynchronize();
}

//====================================================================
void copy_from_qdw(real_t* std_v, const real_t* qdw_v, int nvol, int nin4)
{
  real_t* std_dev = (real_t*)dev_ptr(std_v);
  real_t* qdw_dev = (real_t*)dev_ptr(const_cast<real_t*>(qdw_v));

  int nth = VECTOR_LENGTH;
  int nbl = CEIL_NWP(nvol) / nth;
  if (nbl == 0) nbl = 1;

  copy_from_qdw_kernel<<<nbl, nth>>>(std_dev, qdw_dev, nin4, nvol);
  cudaDeviceSynchronize();
}

//====================================================================
__global__ void copy_kernel(real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in = 0; in < nin; ++in)
  {
    v[IDX2(nin, in, site)] = w[IDX2(nin, in, site)];
  }
}

//====================================================================
void copy(real_t *y, real_t *x, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  copy_kernel<<<nbl, nth>>>(y_dev, x_dev, nin);

  cudaDeviceSynchronize();
}

//====================================================================
void copy(real_t *y, int nv1, real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  copy_kernel<<<nbl, nth>>>(&y_dev[nv1], &x_dev[nv2], nin);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void scal_kernel(real_t *v, real_t a, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in = 0; in < nin; ++in)
  {
    v[IDX2(nin, in, site)] *= a;
  }
}

//====================================================================
void scal(real_t *y, int nv1, real_t a, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  scal_kernel<<<nbl, nth>>>(&y_dev[nv1], a, nin);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void scal_kernel(real_t *v, real_t ar, real_t ai, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int kr = 2 * in2;
    int ki = 2 * in2 + 1;
    real_t vr = v[IDX2(nin, kr, site)];
    real_t vi = v[IDX2(nin, ki, site)];
    v[IDX2(nin, kr, site)] = ar * vr - ai * vi;
    v[IDX2(nin, ki, site)] = ai * vr + ar * vi;
  }
}

//====================================================================
void scal(real_t *y, int nv1, real_t ar, real_t ai, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  scal_kernel<<<nbl, nth>>>(&y_dev[nv1], ar, ai, nin);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void axpy_kernel(real_t *v, real_t a, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in = 0; in < nin; ++in)
  {
    v[IDX2(nin, in, site)] += a * w[IDX2(nin, in, site)];
  }
}

//====================================================================
__global__ void axpy_kernel_opt(real_t *v, real_t a, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int TILE_IN = NWP;
  extern __shared__ real_t buf[];

  real_t *tile = buf + threadIdx.x * TILE_IN;

  for (int in_base = 0; in_base < nin; in_base += TILE_IN)
  {
    int this_tile = min(TILE_IN, nin - in_base);

    for (int j = 0; j < this_tile; ++j)
    {
      int in = in_base + j;
      tile[j] = w[IDX2(nin, in, site)];
    }

    __syncthreads();

    for (int j = 0; j < this_tile; ++j)
    {
      int in = in_base + j;
      v[IDX2(nin, in, site)] += a * tile[j];
    }

    __syncthreads();
  }
}

//====================================================================
// QDW-enhanced axpy kernel: uses dw math for better precision (v = a*w + v)
__global__ void axpy_kernel_qdw(real_t *v, real_t a, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    real_t wh_re = w[base + 0];
    real_t wh_im = w[base + 1];
    real_t wl_re = w[base + 2];
    real_t wl_im = w[base + 3];

    real_t vh_re = v[base + 0];
    real_t vh_im = v[base + 1];
    real_t vl_re = v[base + 2];
    real_t vl_im = v[base + 3];

    real_t th_re, tl_re;
    dw_scal(a, wh_re, wl_re, th_re, tl_re);
    real_t th_im, tl_im;
    dw_scal(a, wh_im, wl_im, th_im, tl_im);

    real_t out_h_re, out_l_re;
    dw_add(vh_re, vl_re, th_re, tl_re, out_h_re, out_l_re);
    real_t out_h_im, out_l_im;
    dw_add(vh_im, vl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = out_h_re;
    v[base + 1] = out_h_im;
    v[base + 2] = out_l_re;
    v[base + 3] = out_l_im;
  }
}

//====================================================================
// DD-pair axpy: v = a*w + v  where a = (a_h, a_l) is a float-pair scalar.
// Uses dw_mul on the scalar-pair × vector-pair so all 48-bit precision in
// the scalar is preserved through the multiplication.
__global__ void axpy_kernel_qdw_pair(real_t *v, real_t a_h, real_t a_l, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    real_t wh_re = w[base + 0];
    real_t wh_im = w[base + 1];
    real_t wl_re = w[base + 2];
    real_t wl_im = w[base + 3];

    real_t vh_re = v[base + 0];
    real_t vh_im = v[base + 1];
    real_t vl_re = v[base + 2];
    real_t vl_im = v[base + 3];

    // (a_h,a_l) * (w_hi,w_lo) — full DD multiply
    real_t th_re, tl_re;
    dw_mul(a_h, a_l, wh_re, wl_re, th_re, tl_re);
    real_t th_im, tl_im;
    dw_mul(a_h, a_l, wh_im, wl_im, th_im, tl_im);

    real_t out_h_re, out_l_re;
    dw_add(vh_re, vl_re, th_re, tl_re, out_h_re, out_l_re);
    real_t out_h_im, out_l_im;
    dw_add(vh_im, vl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = out_h_re;
    v[base + 1] = out_h_im;
    v[base + 2] = out_l_re;
    v[base + 3] = out_l_im;
  }
}

void axpy_pair(real_t *y, int nv1, real_t a_h, real_t a_l,
               real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    axpy_kernel_qdw_pair<<<nbl, nth>>>(&y_dev[nv1], a_h, a_l, &x_dev[nv2], nin);
  } else {
    // Non-DW path collapses the pair into a single float.
    real_t a = a_h + a_l;
    size_t sharedMemSize = nth * NWP * sizeof(real_t);
    axpy_kernel_opt<<<nbl, nth, sharedMemSize>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}

//====================================================================
// Forward declarations: the TW (mode==2) kernels are defined later in this TU.
__global__ void axpy_kernel_tw3(real_t *v, real_t a, real_t *w, int nin);
__global__ void aypx_kernel_tw3(real_t a, real_t *v, real_t *w, int nin);

void axpy(real_t *y, int nv1, real_t a,
          real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    axpy_kernel_qdw<<<nbl, nth>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  } else if (mode == 2) {
    // TW: genuine triple-word a*x + y (was silently single-float before).
    axpy_kernel_tw3<<<nbl, nth>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  } else {
    size_t sharedMemSize = nth * NWP * sizeof(real_t);
    axpy_kernel_opt<<<nbl, nth, sharedMemSize>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}


//====================================================================
__global__ void axpy_kernel(real_t *v, real_t ar, real_t ai,
                            real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int kr = 2 * in2;
    int ki = 2 * in2 + 1;
    real_t wr = w[IDX2(nin, kr, site)];
    real_t wi = w[IDX2(nin, ki, site)];
    v[IDX2(nin, kr, site)] += ar * wr - ai * wi;
    v[IDX2(nin, ki, site)] += ai * wr + ar * wi;
  }
}

//====================================================================
// QDW-enhanced complex axpy kernel: uses __fma_rn
__global__ void axpy_kernel_qdw(real_t *v, real_t ar, real_t ai,
                                real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int kr = 2 * in2;
    int ki = 2 * in2 + 1;
    real_t wr = w[IDX2(nin, kr, site)];
    real_t wi = w[IDX2(nin, ki, site)];
    // v_r += ar*wr - ai*wi, v_i += ai*wr + ar*wi
    v[IDX2(nin, kr, site)] = __fma_rn(ar, wr, __fma_rn(-ai, wi, v[IDX2(nin, kr, site)]));
    v[IDX2(nin, ki, site)] = __fma_rn(ai, wr, __fma_rn( ar, wi, v[IDX2(nin, ki, site)]));
  }
}

//====================================================================
void axpy(real_t *y, int nv1, real_t ar, real_t ai,
          real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    axpy_kernel_qdw<<<nbl, nth>>>(&y_dev[nv1], ar, ai, &x_dev[nv2], nin);
  } else {
    axpy_kernel<<<nbl, nth>>>(&y_dev[nv1], ar, ai, &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void aypx_kernel(real_t a, real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in = 0; in < nin; ++in)
  {
    v[IDX2(nin, in, site)] = a * v[IDX2(nin, in, site)] + w[IDX2(nin, in, site)];
  }
}

//====================================================================
// QDW-enhanced aypx kernel: uses dw math for better precision (v = a*v + w)
__global__ void aypx_kernel_qdw(real_t a, real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    real_t vh_re = v[base + 0];
    real_t vh_im = v[base + 1];
    real_t vl_re = v[base + 2];
    real_t vl_im = v[base + 3];

    real_t wh_re = w[base + 0];
    real_t wh_im = w[base + 1];
    real_t wl_re = w[base + 2];
    real_t wl_im = w[base + 3];

    real_t th_re, tl_re;
    dw_scal(a, vh_re, vl_re, th_re, tl_re);
    real_t th_im, tl_im;
    dw_scal(a, vh_im, vl_im, th_im, tl_im);

    real_t out_h_re, out_l_re;
    dw_add(wh_re, wl_re, th_re, tl_re, out_h_re, out_l_re);
    real_t out_h_im, out_l_im;
    dw_add(wh_im, wl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = out_h_re;
    v[base + 1] = out_h_im;
    v[base + 2] = out_l_re;
    v[base + 3] = out_l_im;
  }
}

//====================================================================
// DD-pair aypx: v = a*v + w, a = (a_h, a_l) float-pair scalar.
__global__ void aypx_kernel_qdw_pair(real_t a_h, real_t a_l, real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    real_t vh_re = v[base + 0];
    real_t vh_im = v[base + 1];
    real_t vl_re = v[base + 2];
    real_t vl_im = v[base + 3];

    real_t wh_re = w[base + 0];
    real_t wh_im = w[base + 1];
    real_t wl_re = w[base + 2];
    real_t wl_im = w[base + 3];

    real_t th_re, tl_re;
    dw_mul(a_h, a_l, vh_re, vl_re, th_re, tl_re);
    real_t th_im, tl_im;
    dw_mul(a_h, a_l, vh_im, vl_im, th_im, tl_im);

    real_t out_h_re, out_l_re;
    dw_add(wh_re, wl_re, th_re, tl_re, out_h_re, out_l_re);
    real_t out_h_im, out_l_im;
    dw_add(wh_im, wl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = out_h_re;
    v[base + 1] = out_h_im;
    v[base + 2] = out_l_re;
    v[base + 3] = out_l_im;
  }
}

void aypx_pair(real_t a_h, real_t a_l, real_t *y, int nv1,
               real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    aypx_kernel_qdw_pair<<<nbl, nth>>>(a_h, a_l, &y_dev[nv1], &x_dev[nv2], nin);
  } else {
    real_t a = a_h + a_l;
    aypx_kernel<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}

//====================================================================
void aypx(real_t a, real_t *y, int nv1,
          real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    aypx_kernel_qdw<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);
  } else if (mode == 2) {
    // TW: genuine triple-word a*y + x (was silently single-float before).
    aypx_kernel_tw3<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);
  } else {
    aypx_kernel<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void aypx_kernel(real_t ar, real_t ai, real_t *v,
                            real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int kr = 2 * in2;
    int ki = 2 * in2 + 1;
    real_t vr = v[IDX2(nin, kr, site)];
    real_t vi = v[IDX2(nin, ki, site)];
    v[IDX2(nin, kr, site)] = ar * vr - ai * vi + w[IDX2(nin, kr, site)];
    v[IDX2(nin, ki, site)] = ai * vr + ar * vi + w[IDX2(nin, ki, site)];
  }
}

//====================================================================
// QDW-enhanced complex aypx kernel: uses dw math for better precision (v = a*v + w)
// a = ar + i ai
// a*v = (ar*vr - ai*vi) + i(ai*vr + ar*vi)
__global__ void aypx_kernel_qdw(real_t ar, real_t ai, real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;

  for (int in4 = 0; in4 < nin / 4; ++in4)
  {
    real_t vh_re = v[IDX2(nin, 4 * in4 + 0, site)];
    real_t vh_im = v[IDX2(nin, 4 * in4 + 1, site)];
    real_t vl_re = v[IDX2(nin, 4 * in4 + 2, site)];
    real_t vl_im = v[IDX2(nin, 4 * in4 + 3, site)];

    real_t wh_re = w[IDX2(nin, 4 * in4 + 0, site)];
    real_t wh_im = w[IDX2(nin, 4 * in4 + 1, site)];
    real_t wl_re = w[IDX2(nin, 4 * in4 + 2, site)];
    real_t wl_im = w[IDX2(nin, 4 * in4 + 3, site)];

    // real part: ar*vr - ai*vi
    real_t p1_r, e1_r; dw_scal(ar, vh_re, vl_re, p1_r, e1_r);
    real_t p2_r, e2_r; dw_scal(ai, vh_im, vl_im, p2_r, e2_r);
    real_t r_h, r_l;   dw_add(p1_r, e1_r, -p2_r, -e2_r, r_h, r_l);

    // imag part: ai*vr + ar*vi
    real_t p1_i, e1_i; dw_scal(ai, vh_re, vl_re, p1_i, e1_i);
    real_t p2_i, e2_i; dw_scal(ar, vh_im, vl_im, p2_i, e2_i);
    real_t i_h, i_l;   dw_add(p1_i, e1_i, p2_i, e2_i, i_h, i_l);

    // add w: v = a*v + w
    real_t out_h_re, out_l_re; dw_add(wh_re, wl_re, r_h, r_l, out_h_re, out_l_re);
    real_t out_h_im, out_l_im; dw_add(wh_im, wl_im, i_h, i_l, out_h_im, out_l_im);

    v[IDX2(nin, 4 * in4 + 0, site)] = out_h_re;
    v[IDX2(nin, 4 * in4 + 1, site)] = out_h_im;
    v[IDX2(nin, 4 * in4 + 2, site)] = out_l_re;
    v[IDX2(nin, 4 * in4 + 3, site)] = out_l_im;
  }
}

//====================================================================
void aypx(real_t ar, real_t ai, real_t *y, int nv1,
          real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    aypx_kernel_qdw<<<nbl, nth>>>(ar, ai, &y_dev[nv1], &x_dev[nv2], nin);
  } else {
    aypx_kernel<<<nbl, nth>>>(ar, ai, &y_dev[nv1], &x_dev[nv2], nin);
  }

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void xI_kernel(real_t *v, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int kr = 2 * in2;
    int ki = 2 * in2 + 1;
    real_t vr = v[IDX2(nin, kr, site)];
    real_t vi = v[IDX2(nin, ki, site)];
    v[IDX2(nin, kr, site)] = -vi;
    v[IDX2(nin, ki, site)] = vr;
  }
}

//====================================================================
void xI(real_t *x, int nin, int nst)
{
  int nth = VECTOR_LENGTH;
  int nbl = nst / nth;

  real_t *x_dev = (real_t *)dev_ptr(x);

  xI_kernel<<<nbl, nth>>>(x_dev, nin);

  cudaDeviceSynchronize();
}

//====================================================================
__global__ void conjg_kernel(real_t *v, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  for (int in2 = 0; in2 < nin2; ++in2)
  {
    int ki = 2 * in2 + 1;
    real_t vi = v[IDX2(nin, ki, site)];
    v[IDX2(nin, ki, site)] = -vi;
  }
}

//====================================================================
void conjg(real_t *x, int nin, int nst)
{
  int nth = VECTOR_LENGTH;
  int nbl = nst / nth;

  real_t *x_dev = (real_t *)dev_ptr(x);

  conjg_kernel<<<nbl, nth>>>(x_dev, nin);

  cudaDeviceSynchronize();
}

//-------------------
__global__ void reduce_kernel_multiblocks(real_t *red, int nvol)
{
  extern __shared__ real_t sdata[];

  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum = 0.0;

  for (int i = idx; i < nvol; i += gridSize)
  {
    sum += red[i];
  }

  sdata[tid] = sum;
  __syncthreads();

  for (int s = blockDim.x / 2; s > 0; s >>= 1)
  {
    if (tid < s)
    {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  if (tid == 0)
  {
    red[blockIdx.x] = sdata[0];
  }
}

//====================================================================
__global__ void reduce_kernel(real_t *red, int nvol)
{
  __shared__ real_t red2[VECTOR_LENGTH];

  const int ith = threadIdx.x;
  const int nth = VECTOR_LENGTH;
  const int nvol2 = nvol / nth;

  real_t at = 0.0;
  for (int ist2 = 0; ist2 < nvol2; ++ist2)
  {
    int ist = ith + nth * ist2;
    at += red[ist];
  }
  red2[ith] = at;

  __syncthreads();

  if (ith == 0)
  {
    real_t at2 = red2[0];
#pragma unroll
    for (int i = 1; i < VECTOR_LENGTH; ++i)
    {
      at2 += red2[i];
    }
    red[0] = at2;
  }
}

//====================================================================
__global__ void norm2_kernel(real_t *red, real_t *v1,
                             int nin, int nvol, int nex)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    real_t at = 0.0;
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = ist + nvol * ex;
#pragma unroll
      for (int in = 0; in < nin; ++in)
      {
        real_t vt1 = v1[IDX2(nin, in, ist2)];
        at += vt1 * vt1;
      }
    }
    red[ist] = at;
  }
}
//====================================================================
__global__ void norm2_reduce_fused_kernel(real_t *red,
                                          const real_t *__restrict__ v1,
                                          int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];

  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum = 0.0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    real_t at = 0.0;
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      for (int in = 0; in < nin; ++in)
      {
        real_t vt = v1[IDX2(nin, in, ist2)];
        at += vt * vt;
      }
    }
    sum += at;
  }

  sdata[tid] = sum;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {

    sum = sdata[tid];

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      sum += __shfl_down_sync(0xffffffff, sum, offset);
    }

    if (tid == 0)
    {
      red[blockIdx.x] = sum;
    }
  }
}

//====================================================================
// Generic (real_t-typed) TwoSum/TwoProd for QDW and QTW kernels.
// For real_t=double: exact error-free arithmetic.
// For real_t=float:  F32 via __fmaf_rn, giving double-float precision at FP32 throughput.
__device__ __forceinline__ static void TwoSum_r(real_t a, real_t b, real_t &s, real_t &e) {
    s = a + b;
    real_t v = s - a;
    e = (a - (s - v)) + (b - v);
}
__device__ __forceinline__ static void TwoProd_r(real_t a, real_t b, real_t &p, real_t &e) {
    p = a * b;
    e = fma(a, b, -p);
}

//====================================================================
// QDW-enhanced norm2: Kahan compensated summation in thread-level accumulation
__global__ void norm2_reduce_fused_kernel_qdw(real_t *red,
                                              const real_t *__restrict__ v1,
                                              int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];

  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum  = 0.0;
  real_t comp = 0.0;  // Kahan compensation

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      
      const int nin4 = nin / 4;
      for (int in4 = 0; in4 < nin4; ++in4)
      {
        int base = 4 * IDX2(nin4, in4, ist2);
        real_t h_re = v1[base + 0];
        real_t h_im = v1[base + 1];
        real_t l_re = v1[base + 2];
        real_t l_im = v1[base + 3];

        // Exact real square
        real_t p1, e1;
        TwoProd(h_re, h_re, p1, e1);
        real_t lo = real_t(2.0) * h_re * l_re + e1;

        real_t dh, dl;
        TwoSum(sum, p1, dh, dl);
        real_t temp_sum  = dh;
        real_t temp_comp = comp + dl + lo;
        real_t sum_out, comp_out;
        TwoSum(temp_sum, temp_comp, sum_out, comp_out);
        sum  = sum_out;
        comp = comp_out;

        // Exact imaginary square
        TwoProd(h_im, h_im, p1, e1);
        lo = real_t(2.0) * h_im * l_im + e1;

        TwoSum(sum, p1, dh, dl);
        temp_sum  = dh;
        temp_comp = comp + dl + lo;
        TwoSum(temp_sum, temp_comp, sum_out, comp_out);
        sum  = sum_out;
        comp = comp_out;
      }
    }
  }

  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;
  real_t nh, nl;
  TwoSum_r(sum, comp, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_hi[tid], al = sh_lo[tid];
    real_t bh, bl;

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      bh = __shfl_down_sync(0xffffffff, ah, offset);
      bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      ah = rh;
      al = al + bl + t;
    }

    if (tid == 0)
    {
      red[blockIdx.x]                        = ah;
      red[blockIdx.x + MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
// QTW-enhanced norm2: dual-word (TwoSum_r/TwoProd_r) thread accumulator.
// Works in real_t precision: F32-QTW for float, FP64-QTW for double.
// Outputs high word to red[blockIdx.x], low word to red[blockIdx.x + MAX_THREAD_PER_BLOCK].
__global__ void norm2_reduce_fused_kernel_qtw(real_t *red,
                                              const real_t *__restrict__ v1,
                                              int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;

  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum_h = 0, sum_l = 0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      for (int in = 0; in < nin; ++in)
      {
        real_t vt = v1[IDX2(nin, in, ist2)];
        real_t p, e;
        TwoProd_r(vt, vt, p, e);
        real_t sh, st;
        TwoSum_r(sum_h, p, sh, st);
        sum_h = sh;
        sum_l = sum_l + e + st;
      }
    }
  }

  // Normalize dual-word before storing to shared memory
  real_t nh, nl;
  TwoSum_r(sum_h, sum_l, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_hi[tid], al = sh_lo[tid];
    real_t bh, bl;

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      bh = __shfl_down_sync(0xffffffff, ah, offset);
      bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      ah = rh;
      al = al + bl + t;
    }

    if (tid == 0)
    {
      red[blockIdx.x]                        = ah;
      red[blockIdx.x + MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
// Second-pass reducer for QTW dual-word partial sums.
// hi part: red[0..n-1], lo part: red[MAX_THREAD_PER_BLOCK..MAX_THREAD_PER_BLOCK+n-1].
__global__ void reduce_kernel_multiblocks_qtw(real_t *red, int n)
{
  extern __shared__ real_t sdata[];
  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;

  const int tid = threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum_h = 0, sum_l = 0;

  for (int i = tid; i < n; i += gridSize)
  {
    real_t ah = red[i];
    real_t al = red[i + MAX_THREAD_PER_BLOCK];
    real_t rh, t;
    TwoSum_r(sum_h, ah, rh, t);
    sum_h = rh;
    sum_l = sum_l + al + t;
  }

  real_t nh, nl;
  TwoSum_r(sum_h, sum_l, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid == 0)
  {
    red[0] = sh_hi[0] + sh_lo[0];
  }
}

//====================================================================
// QTW/QDW final-pass reducer that preserves the dual-word pair.
// Writes hi to red[0] and lo to red[1] (instead of collapsing). Used by the
// _pair host wrappers so the CG layer can carry DD scalars without collapsing
// to single real_t at the device→host boundary.
__global__ void reduce_kernel_multiblocks_qtw_pair(real_t *red, int n)
{
  extern __shared__ real_t sdata[];
  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;

  const int tid = threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum_h = 0, sum_l = 0;

  for (int i = tid; i < n; i += gridSize)
  {
    real_t ah = red[i];
    real_t al = red[i + MAX_THREAD_PER_BLOCK];
    real_t rh, t;
    TwoSum_r(sum_h, ah, rh, t);
    sum_h = rh;
    sum_l = sum_l + al + t;
  }

  real_t nh, nl;
  TwoSum_r(sum_h, sum_l, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid == 0)
  {
    real_t nh2, nl2;
    TwoSum_r(sh_hi[0], sh_lo[0], nh2, nl2);
    red[0] = nh2;
    red[1] = nl2;
  }
}

//====================================================================
real_t norm2(real_t *x1, real_t *red1, int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  int nth2 = VECTOR_LENGTH;
  int nbl2 = 1;

  int threadsPerBlock = nth2;
  int blocksPerGrid = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    norm2_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else if (mode == 2) {
    norm2_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else {
    norm2_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t a;
  CHECK(cudaMemcpy(&a, &red1_dev[0], sizeof(real_t),
                   cudaMemcpyDeviceToHost));

  return a;
}

//====================================================================
// DD-pair norm2: returns (h, l) preserving the dual-word accumulator. Only
// meaningful for mode 1/2; for mode 0 we fall back to the single-real_t path
// and return l=0.
void norm2_pair(real_t *h, real_t *l, real_t *x1, real_t *red1,
                int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    norm2_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else if (mode == 2) {
    norm2_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else {
    norm2_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
    reduce_kernel_multiblocks<<<1, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t pair[2] = {0, 0};
  if (mode == 1 || mode == 2) {
    CHECK(cudaMemcpy(pair, &red1_dev[0], 2 * sizeof(real_t), cudaMemcpyDeviceToHost));
  } else {
    CHECK(cudaMemcpy(&pair[0], &red1_dev[0], sizeof(real_t), cudaMemcpyDeviceToHost));
  }
  *h = pair[0];
  *l = pair[1];
}

//====================================================================
__global__ void dot_kernel(real_t *red, real_t *v1, real_t *v2,
                           int nin, int nvol, int nex)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;

  real_t at = 0.0;
  for (int ex = 0; ex < nex; ++ex)
  {
    int ist2 = ist + nvol * ex;
    for (int in = 0; in < nin; ++in)
    {
      at += v1[IDX2(nin, in, ist2)] * v2[IDX2(nin, in, ist2)];
    }
  }

  red[ist] = at;
}
//===================================================================
__global__ void dot_reduce_fused_kernel(real_t *red,
                                        const real_t *__restrict__ v1,
                                        const real_t *__restrict__ v2,
                                        int nin, int nvol, int nex)
{

  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  extern __shared__ real_t sdata[];

  real_t at = 0.0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      for (int in = 0; in < nin; ++in)
      {
        at += v1[IDX2(nin, in, ist2)] * v2[IDX2(nin, in, ist2)];
      }
    }
  }

  sdata[tid] = at;
  __syncthreads();

  for (int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {

    at = sdata[tid];

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      at += __shfl_down_sync(0xffffffff, at, offset);
    }

    if (tid == 0)
    {
      red[blockIdx.x] = at;
    }
  }
}

//====================================================================
// QDW-enhanced dot: Kahan compensated summation
__global__ void dot_reduce_fused_kernel_qdw(real_t *red,
                                            const real_t *__restrict__ v1,
                                            const real_t *__restrict__ v2,
                                            int nin, int nvol, int nex)
{
  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  extern __shared__ real_t sdata[];

  real_t at   = 0.0;
  real_t comp = 0.0;  // Kahan compensation

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      const int nin4 = nin / 4;
      for (int in4 = 0; in4 < nin4; ++in4)
      {
        int base = 4 * IDX2(nin4, in4, ist2);
        real_t v1h_re = v1[base + 0];
        real_t v1h_im = v1[base + 1];
        real_t v1l_re = v1[base + 2];
        real_t v1l_im = v1[base + 3];

        real_t v2h_re = v2[base + 0];
        real_t v2h_im = v2[base + 1];
        real_t v2l_re = v2[base + 2];
        real_t v2l_im = v2[base + 3];

        // Real dot product part: v1.re * v2.re + v1.im * v2.im
        real_t p1_r, e1_r, p2_r, e2_r;
        TwoProd(v1h_re, v2h_re, p1_r, e1_r);
        TwoProd(v1h_im, v2h_im, p2_r, e2_r);

        real_t dh_r, dl_r;
        TwoSum(p1_r, p2_r, dh_r, dl_r);

        real_t lo_r = v1h_re * v2l_re + v1l_re * v2h_re + v1h_im * v2l_im + v1l_im * v2h_im + e1_r + e2_r + dl_r;

        real_t dh, dl;
        TwoSum(at, dh_r, dh, dl);
        real_t temp_sum  = dh;
        real_t temp_comp = comp + dl + lo_r;
        real_t at_out, comp_out;
        TwoSum(temp_sum, temp_comp, at_out, comp_out);
        at   = at_out;
        comp = comp_out;
      }
    }
  }

  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;
  real_t nh, nl;
  TwoSum_r(at, comp, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_hi[tid], al = sh_lo[tid];
    real_t bh, bl;

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      bh = __shfl_down_sync(0xffffffff, ah, offset);
      bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      ah = rh;
      al = al + bl + t;
    }

    if (tid == 0)
    {
      red[blockIdx.x]                        = ah;
      red[blockIdx.x + MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
// QTW-enhanced dot: dual-word thread accumulator for standard double fields.
__global__ void dot_reduce_fused_kernel_qtw(real_t *red,
                                            const real_t *__restrict__ v1,
                                            const real_t *__restrict__ v2,
                                            int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  real_t *sh_hi = sdata;
  real_t *sh_lo = sdata + blockDim.x;

  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sum_h = 0, sum_l = 0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      for (int in = 0; in < nin; ++in)
      {
        real_t a = v1[IDX2(nin, in, ist2)];
        real_t b = v2[IDX2(nin, in, ist2)];
        real_t p, e;
        TwoProd_r(a, b, p, e);
        real_t rh, t;
        TwoSum_r(sum_h, p, rh, t);
        sum_h = rh;
        sum_l = sum_l + e + t;
      }
    }
  }

  real_t nh, nl;
  TwoSum_r(sum_h, sum_l, nh, nl);
  sh_hi[tid] = nh;
  sh_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah = sh_hi[tid],   al = sh_lo[tid];
      real_t bh = sh_hi[tid+s], bl = sh_lo[tid+s];
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      sh_hi[tid] = rh;
      sh_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_hi[tid], al = sh_lo[tid];
    real_t bh, bl;
#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      bh = __shfl_down_sync(0xffffffff, ah, offset);
      bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, t;
      TwoSum_r(ah, bh, rh, t);
      ah = rh;
      al = al + bl + t;
    }
    if (tid == 0)
    {
      red[blockIdx.x]                        = ah;
      red[blockIdx.x + MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
real_t dot(real_t *x1, real_t *x2, real_t *red1,
           int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *x2_dev = (real_t *)dev_ptr(x2);

  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  int nth2 = VECTOR_LENGTH;
  int nbl2 = 1;
  int nSM = getNumSMs();
  int threadsPerBlock = nth2;
  int blocksPerGrid = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    dot_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev,
                                                                                     x1_dev, x2_dev,
                                                                                     nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else if (mode == 2) {
    dot_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev,
                                                                                     x1_dev, x2_dev,
                                                                                     nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else {
    dot_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev,
                                                                              x1_dev, x2_dev,
                                                                              nin, nvol, nex);
    reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t a;
  CHECK(cudaMemcpy(&a, &red1_dev[0], sizeof(real_t),
                   cudaMemcpyDeviceToHost));

  return a;
}

//====================================================================
// DD-pair dot: returns (h, l) preserving the dual-word accumulator.
void dot_pair(real_t *h, real_t *l, real_t *x1, real_t *x2, real_t *red1,
              int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *x2_dev   = (real_t *)dev_ptr(x2);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    dot_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else if (mode == 2) {
    dot_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
  } else {
    dot_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks<<<1, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t pair[2] = {0, 0};
  if (mode == 1 || mode == 2) {
    CHECK(cudaMemcpy(pair, &red1_dev[0], 2 * sizeof(real_t), cudaMemcpyDeviceToHost));
  } else {
    CHECK(cudaMemcpy(&pair[0], &red1_dev[0], sizeof(real_t), cudaMemcpyDeviceToHost));
  }
  *h = pair[0];
  *l = pair[1];
}

//====================================================================
__global__ void dotc_kernel(real_t *red1, real_t *red2,
                            real_t *v1, real_t *v2,
                            int nin, int nvol, int nex)
{
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin2 = nin / 2;

  real_t ar = 0.0;
  real_t ai = 0.0;

  for (int ex = 0; ex < nex; ++ex)
  {
    int ist2 = ist + nvol * ex;

    for (int in2 = 0; in2 < nin2; ++in2)
    {
      int kr = 2 * in2;
      int ki = 2 * in2 + 1;

      real_t v1r = v1[IDX2(nin, kr, ist2)];
      real_t v1i = v1[IDX2(nin, ki, ist2)];
      real_t v2r = v2[IDX2(nin, kr, ist2)];
      real_t v2i = v2[IDX2(nin, ki, ist2)];

      ar += v1r * v2r + v1i * v2i;
      ai += v1r * v2i - v1i * v2r;
    }
  }

  red1[ist] = ar;
  red2[ist] = ai;
}
//===================================================================
__global__ void dotc_reduce_fused_kernel(real_t *red1, real_t *red2,
                                         const real_t *__restrict__ v1,
                                         const real_t *__restrict__ v2,
                                         int nin, int nvol, int nex)
{

  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  extern __shared__ real_t sdata[];
  real_t *sh_r = sdata;
  real_t *sh_i = sdata + blockDim.x;

  const int nin2 = nin / 2;

  real_t ar = 0.0;
  real_t ai = 0.0;

  for (int ist = idx; ist < nvol; ist += gridSize)
  {

    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = ist + nvol * ex;

      for (int in2 = 0; in2 < nin2; ++in2)
      {
        int kr = 2 * in2;
        int ki = 2 * in2 + 1;

        real_t v1r = v1[IDX2(nin, kr, ist2)];
        real_t v1i = v1[IDX2(nin, ki, ist2)];
        real_t v2r = v2[IDX2(nin, kr, ist2)];
        real_t v2i = v2[IDX2(nin, ki, ist2)];

        ar += v1r * v2r + v1i * v2i;
        ai += v1r * v2i - v1i * v2r;
      }
    }
  }

  sh_r[tid] = ar;
  sh_i[tid] = ai;
  __syncthreads();

  for (int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      sh_r[tid] += sh_r[tid + s];
      sh_i[tid] += sh_i[tid + s];
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {

    ar = sh_r[tid];
    ai = sh_i[tid];

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      ar += __shfl_down_sync(0xffffffff, ar, offset);
      ai += __shfl_down_sync(0xffffffff, ai, offset);
    }

    if (tid == 0)
    {
      red1[blockIdx.x] = ar;
      red2[blockIdx.x] = ai;
    }
  }
}
//====================================================================
//====================================================================
// QDW-enhanced dotc: Kahan compensated summation for complex dot product
__global__ void dotc_reduce_fused_kernel_qdw(real_t *red1, real_t *red2,
                                             const real_t *__restrict__ v1,
                                             const real_t *__restrict__ v2,
                                             int nin, int nvol, int nex)
{
  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  extern __shared__ real_t sdata[];
  real_t *sh_r_hi = sdata;
  real_t *sh_r_lo = sdata +     blockDim.x;
  real_t *sh_i_hi = sdata + 2 * blockDim.x;
  real_t *sh_i_lo = sdata + 3 * blockDim.x;

  const int nin2 = nin / 2;

  real_t ar = 0.0, ai = 0.0;
  real_t comp_r = 0.0, comp_i = 0.0;

  for (int ist = idx; ist < nvol; ist += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = ist + nvol * ex;

      const int nin4 = nin / 4;
      for (int in4 = 0; in4 < nin4; ++in4)
      {
        int base = 4 * IDX2(nin4, in4, ist2);
        real_t v1h_re = v1[base + 0];
        real_t v1h_im = v1[base + 1];
        real_t v1l_re = v1[base + 2];
        real_t v1l_im = v1[base + 3];

        real_t v2h_re = v2[base + 0];
        real_t v2h_im = v2[base + 1];
        real_t v2l_re = v2[base + 2];
        real_t v2l_im = v2[base + 3];

        // Real part: v1.re * v2.re + v1.im * v2.im
        real_t p1_r, e1_r, p2_r, e2_r;
        TwoProd(v1h_re, v2h_re, p1_r, e1_r);
        TwoProd(v1h_im, v2h_im, p2_r, e2_r);

        real_t dh_r, dl_r;
        TwoSum(p1_r, p2_r, dh_r, dl_r);

        real_t lo_r = v1h_re * v2l_re + v1l_re * v2h_re + v1h_im * v2l_im + v1l_im * v2h_im + e1_r + e2_r + dl_r;

        real_t dh, dl;
        TwoSum(ar, dh_r, dh, dl);
        real_t temp_sum_r  = dh;
        real_t temp_comp_r = comp_r + dl + lo_r;
        real_t ar_out, comp_r_out;
        TwoSum(temp_sum_r, temp_comp_r, ar_out, comp_r_out);
        ar     = ar_out;
        comp_r = comp_r_out;

        // Imag part: v1.re * v2.im - v1.im * v2.re
        real_t p1_i, e1_i, p2_i, e2_i;
        TwoProd(v1h_re, v2h_im, p1_i, e1_i);
        TwoProd(-v1h_im, v2h_re, p2_i, e2_i);

        real_t dh_i, dl_i;
        TwoSum(p1_i, p2_i, dh_i, dl_i);

        real_t lo_i = v1h_re * v2l_im + v1l_re * v2h_im - v1h_im * v2l_re - v1l_im * v2h_re + e1_i + e2_i + dl_i;

        TwoSum(ai, dh_i, dh, dl);
        real_t temp_sum_i  = dh;
        real_t temp_comp_i = comp_i + dl + lo_i;
        real_t ai_out, comp_i_out;
        TwoSum(temp_sum_i, temp_comp_i, ai_out, comp_i_out);
        ai     = ai_out;
        comp_i = comp_i_out;
      }
    }
  }

  real_t nr_h, nr_l, ni_h, ni_l;
  TwoSum_r(ar, comp_r, nr_h, nr_l);
  TwoSum_r(ai, comp_i, ni_h, ni_l);
  sh_r_hi[tid] = nr_h;
  sh_r_lo[tid] = nr_l;
  sh_i_hi[tid] = ni_h;
  sh_i_lo[tid] = ni_l;
  __syncthreads();

  for (int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t arh = sh_r_hi[tid], arl = sh_r_lo[tid];
      real_t brh = sh_r_hi[tid+s], brl = sh_r_lo[tid+s];
      real_t aih = sh_i_hi[tid], ail = sh_i_lo[tid];
      real_t bih = sh_i_hi[tid+s], bil = sh_i_lo[tid+s];
      real_t rh, t;
      TwoSum_r(arh, brh, rh, t);
      sh_r_hi[tid] = rh;
      sh_r_lo[tid] = arl + brl + t;
      TwoSum_r(aih, bih, rh, t);
      sh_i_hi[tid] = rh;
      sh_i_lo[tid] = ail + bil + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t arh = sh_r_hi[tid], arl = sh_r_lo[tid];
    real_t aih = sh_i_hi[tid], ail = sh_i_lo[tid];
    real_t brh, brl, bih, bil;

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      brh = __shfl_down_sync(0xffffffff, arh, offset);
      brl = __shfl_down_sync(0xffffffff, arl, offset);
      bih = __shfl_down_sync(0xffffffff, aih, offset);
      bil = __shfl_down_sync(0xffffffff, ail, offset);
      real_t rh, t;
      TwoSum_r(arh, brh, rh, t);
      arh = rh;
      arl = arl + brl + t;
      TwoSum_r(aih, bih, rh, t);
      aih = rh;
      ail = ail + bil + t;
    }

    if (tid == 0)
    {
      red1[blockIdx.x]                        = arh;
      red1[blockIdx.x + MAX_THREAD_PER_BLOCK] = arl;
      red2[blockIdx.x]                        = aih;
      red2[blockIdx.x + MAX_THREAD_PER_BLOCK] = ail;
    }
  }
}

//====================================================================
// QTW-enhanced dotc: dual-word thread accumulator for complex dot product.
__global__ void dotc_reduce_fused_kernel_qtw(real_t *red1, real_t *red2,
                                             const real_t *__restrict__ v1,
                                             const real_t *__restrict__ v2,
                                             int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  real_t *shr_hi = sdata;
  real_t *shr_lo = sdata +   blockDim.x;
  real_t *shi_hi = sdata + 2*blockDim.x;
  real_t *shi_lo = sdata + 3*blockDim.x;

  const int tid = threadIdx.x;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;
  const int nin2 = nin / 2;

  real_t ar_h = 0, ar_l = 0;
  real_t ai_h = 0, ai_l = 0;

  for (int ist = idx; ist < nvol; ist += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = ist + nvol * ex;
      for (int in2 = 0; in2 < nin2; ++in2)
      {
        int kr = 2 * in2;
        int ki = 2 * in2 + 1;
        real_t v1r = v1[IDX2(nin, kr, ist2)];
        real_t v1i = v1[IDX2(nin, ki, ist2)];
        real_t v2r = v2[IDX2(nin, kr, ist2)];
        real_t v2i = v2[IDX2(nin, ki, ist2)];

        // real part: v1r*v2r + v1i*v2i
        real_t p1, e1, p2, e2;
        TwoProd_r(v1r, v2r, p1, e1);
        TwoProd_r(v1i, v2i, p2, e2);
        real_t ps, pt;
        TwoSum_r(p1, p2, ps, pt);
        real_t rh, rt;
        TwoSum_r(ar_h, ps, rh, rt);
        ar_h = rh;
        ar_l = ar_l + e1 + e2 + pt + rt;

        // imag part: v1r*v2i - v1i*v2r
        real_t p3, e3, p4, e4;
        TwoProd_r(v1r,  v2i, p3, e3);
        TwoProd_r(-v1i, v2r, p4, e4);
        TwoSum_r(p3, p4, ps, pt);
        TwoSum_r(ai_h, ps, rh, rt);
        ai_h = rh;
        ai_l = ai_l + e3 + e4 + pt + rt;
      }
    }
  }

  real_t nh, nl;
  TwoSum_r(ar_h, ar_l, nh, nl);
  shr_hi[tid] = nh;  shr_lo[tid] = nl;
  TwoSum_r(ai_h, ai_l, nh, nl);
  shi_hi[tid] = nh;  shi_lo[tid] = nl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t ah, al, bh, bl, rh, t;
      ah = shr_hi[tid];   al = shr_lo[tid];
      bh = shr_hi[tid+s]; bl = shr_lo[tid+s];
      TwoSum_r(ah, bh, rh, t);
      shr_hi[tid] = rh;  shr_lo[tid] = al + bl + t;

      ah = shi_hi[tid];   al = shi_lo[tid];
      bh = shi_hi[tid+s]; bl = shi_lo[tid+s];
      TwoSum_r(ah, bh, rh, t);
      shi_hi[tid] = rh;  shi_lo[tid] = al + bl + t;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t arh = shr_hi[tid], arl = shr_lo[tid];
    real_t aih = shi_hi[tid], ail = shi_lo[tid];
    real_t bh, bl;

#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      real_t rh, t;
      bh = __shfl_down_sync(0xffffffff, arh, offset);
      bl = __shfl_down_sync(0xffffffff, arl, offset);
      TwoSum_r(arh, bh, rh, t);  arh = rh;  arl = arl + bl + t;

      bh = __shfl_down_sync(0xffffffff, aih, offset);
      bl = __shfl_down_sync(0xffffffff, ail, offset);
      TwoSum_r(aih, bh, rh, t);  aih = rh;  ail = ail + bl + t;
    }

    if (tid == 0)
    {
      red1[blockIdx.x]                        = arh;
      red1[blockIdx.x + MAX_THREAD_PER_BLOCK] = arl;
      red2[blockIdx.x]                        = aih;
      red2[blockIdx.x + MAX_THREAD_PER_BLOCK] = ail;
    }
  }
}

//====================================================================
void dotc(real_t *ar, real_t *ai, real_t *x1, real_t *x2,
          real_t *red1, real_t *red2, int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *x2_dev = (real_t *)dev_ptr(x2);

  real_t *red1_dev = (real_t *)dev_ptr(red1);
  real_t *red2_dev = (real_t *)dev_ptr(red2);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  int nth2 = VECTOR_LENGTH;
  int nbl2 = 1;
  int threadsPerBlock  = nth2;
  int blocksPerGrid    = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    dotc_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 4 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                        x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red2_dev, blocksPerGrid);
  } else if (mode == 2) {
    dotc_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 4 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                        x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks_qtw<<<nbl2, threadsPerBlock, 2*sharedMemSize>>>(red2_dev, blocksPerGrid);
  } else {
    dotc_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, 2 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                    x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red2_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t atr, ati;
  CHECK(cudaMemcpy(&atr, &red1_dev[0], sizeof(real_t),
                   cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(&ati, &red2_dev[0], sizeof(real_t),
                   cudaMemcpyDeviceToHost));

  *ar = atr;
  *ai = ati;
}

//====================================================================
// DD-pair dotc: returns (ar_h, ar_l), (ai_h, ai_l).
void dotc_pair(real_t *ar_h, real_t *ar_l, real_t *ai_h, real_t *ai_l,
               real_t *x1, real_t *x2,
               real_t *red1, real_t *red2,
               int nin, int nvol, int nex, int mode)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *x2_dev   = (real_t *)dev_ptr(x2);
  real_t *red1_dev = (real_t *)dev_ptr(red1);
  real_t *red2_dev = (real_t *)dev_ptr(red2);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  if (mode == 1) {
    dotc_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 4 * sharedMemSize>>>(red1_dev, red2_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red2_dev, blocksPerGrid);
  } else if (mode == 2) {
    dotc_reduce_fused_kernel_qtw<<<blocksPerGrid, threadsPerBlock, 4 * sharedMemSize>>>(red1_dev, red2_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks_qtw_pair<<<1, threadsPerBlock, 2*sharedMemSize>>>(red2_dev, blocksPerGrid);
  } else {
    dotc_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, 2 * sharedMemSize>>>(red1_dev, red2_dev, x1_dev, x2_dev, nin, nvol, nex);
    reduce_kernel_multiblocks<<<1, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
    reduce_kernel_multiblocks<<<1, threadsPerBlock, sharedMemSize>>>(red2_dev, blocksPerGrid);
  }

  cudaDeviceSynchronize();

  real_t r1[2] = {0, 0}, r2[2] = {0, 0};
  if (mode == 1 || mode == 2) {
    CHECK(cudaMemcpy(r1, &red1_dev[0], 2 * sizeof(real_t), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(r2, &red2_dev[0], 2 * sizeof(real_t), cudaMemcpyDeviceToHost));
  } else {
    CHECK(cudaMemcpy(&r1[0], &red1_dev[0], sizeof(real_t), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(&r2[0], &red2_dev[0], sizeof(real_t), cudaMemcpyDeviceToHost));
  }
  *ar_h = r1[0]; *ar_l = r1[1];
  *ai_h = r2[0]; *ai_l = r2[1];
}

//====================================================================
// QTW (triple-word) BLAS kernels. Field layout: 6 real_t per complex,
// {r_hi, i_hi, r_mid, i_mid, r_lo, i_lo}. All arithmetic stays in real_t —
// when real_t = float the chain is FP32-only via ThreeSum/ThreeProd/tw_mul.
//====================================================================

// Helper: read one complex (6 reals starting at base) into 6 named scalars.
#define TW3_LOAD6(prefix, ptr, base) \
    real_t prefix##_rh = (ptr)[(base) + 0]; \
    real_t prefix##_ih = (ptr)[(base) + 1]; \
    real_t prefix##_rm = (ptr)[(base) + 2]; \
    real_t prefix##_im = (ptr)[(base) + 3]; \
    real_t prefix##_rl = (ptr)[(base) + 4]; \
    real_t prefix##_il = (ptr)[(base) + 5];

#define TW3_STORE6(ptr, base, rh, ih, rm, im, rl, il) \
    (ptr)[(base) + 0] = (rh); \
    (ptr)[(base) + 1] = (ih); \
    (ptr)[(base) + 2] = (rm); \
    (ptr)[(base) + 3] = (im); \
    (ptr)[(base) + 4] = (rl); \
    (ptr)[(base) + 5] = (il);

//====================================================================
// axpy_tw3: v = a*w + v   (a = real scalar T)
__global__ void axpy_kernel_tw3(real_t *v, real_t a, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin6 = nin / 6;
  for (int in6 = 0; in6 < nin6; ++in6)
  {
    int base = 6 * IDX2(nin6, in6, site);
    TW3_LOAD6(w, w, base);
    TW3_LOAD6(v, v, base);

    real_t tr_h, tr_m, tr_l;
    tw_scal(a, w_rh, w_rm, w_rl, tr_h, tr_m, tr_l);
    real_t ti_h, ti_m, ti_l;
    tw_scal(a, w_ih, w_im, w_il, ti_h, ti_m, ti_l);

    real_t or_h, or_m, or_l;
    tw_add(v_rh, v_rm, v_rl, tr_h, tr_m, tr_l, or_h, or_m, or_l);
    real_t oi_h, oi_m, oi_l;
    tw_add(v_ih, v_im, v_il, ti_h, ti_m, ti_l, oi_h, oi_m, oi_l);

    TW3_STORE6(v, base, or_h, oi_h, or_m, oi_m, or_l, oi_l);
  }
}

// axpy_tw3 with TW-pair scalar: a = (a_h, a_m, a_l) — full FP32 triple-word
__global__ void axpy_kernel_tw3_pair(real_t *v,
                                     real_t a_h, real_t a_m, real_t a_l,
                                     real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin6 = nin / 6;
  for (int in6 = 0; in6 < nin6; ++in6)
  {
    int base = 6 * IDX2(nin6, in6, site);
    TW3_LOAD6(w, w, base);
    TW3_LOAD6(v, v, base);

    real_t tr_h, tr_m, tr_l;
    tw_mul(a_h, a_m, a_l, w_rh, w_rm, w_rl, tr_h, tr_m, tr_l);
    real_t ti_h, ti_m, ti_l;
    tw_mul(a_h, a_m, a_l, w_ih, w_im, w_il, ti_h, ti_m, ti_l);

    real_t or_h, or_m, or_l;
    tw_add(v_rh, v_rm, v_rl, tr_h, tr_m, tr_l, or_h, or_m, or_l);
    real_t oi_h, oi_m, oi_l;
    tw_add(v_ih, v_im, v_il, ti_h, ti_m, ti_l, oi_h, oi_m, oi_l);

    TW3_STORE6(v, base, or_h, oi_h, or_m, oi_m, or_l, oi_l);
  }
}

void axpy_tw3(real_t *y, int nv1, real_t a,
              real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);
  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  axpy_kernel_tw3<<<nbl, nth>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  cudaDeviceSynchronize();
}

void axpy_tw3_pair(real_t *y, int nv1,
                   real_t a_h, real_t a_m, real_t a_l,
                   real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);
  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  axpy_kernel_tw3_pair<<<nbl, nth>>>(&y_dev[nv1], a_h, a_m, a_l, &x_dev[nv2], nin);
  cudaDeviceSynchronize();
}

//====================================================================
// aypx_tw3: v = a*v + w  (real scalar)
__global__ void aypx_kernel_tw3(real_t a, real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin6 = nin / 6;
  for (int in6 = 0; in6 < nin6; ++in6)
  {
    int base = 6 * IDX2(nin6, in6, site);
    TW3_LOAD6(v, v, base);
    TW3_LOAD6(w, w, base);

    real_t tr_h, tr_m, tr_l;
    tw_scal(a, v_rh, v_rm, v_rl, tr_h, tr_m, tr_l);
    real_t ti_h, ti_m, ti_l;
    tw_scal(a, v_ih, v_im, v_il, ti_h, ti_m, ti_l);

    real_t or_h, or_m, or_l;
    tw_add(w_rh, w_rm, w_rl, tr_h, tr_m, tr_l, or_h, or_m, or_l);
    real_t oi_h, oi_m, oi_l;
    tw_add(w_ih, w_im, w_il, ti_h, ti_m, ti_l, oi_h, oi_m, oi_l);

    TW3_STORE6(v, base, or_h, oi_h, or_m, oi_m, or_l, oi_l);
  }
}

__global__ void aypx_kernel_tw3_pair(real_t a_h, real_t a_m, real_t a_l,
                                     real_t *v, real_t *w, int nin)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  const int nin6 = nin / 6;
  for (int in6 = 0; in6 < nin6; ++in6)
  {
    int base = 6 * IDX2(nin6, in6, site);
    TW3_LOAD6(v, v, base);
    TW3_LOAD6(w, w, base);

    real_t tr_h, tr_m, tr_l;
    tw_mul(a_h, a_m, a_l, v_rh, v_rm, v_rl, tr_h, tr_m, tr_l);
    real_t ti_h, ti_m, ti_l;
    tw_mul(a_h, a_m, a_l, v_ih, v_im, v_il, ti_h, ti_m, ti_l);

    real_t or_h, or_m, or_l;
    tw_add(w_rh, w_rm, w_rl, tr_h, tr_m, tr_l, or_h, or_m, or_l);
    real_t oi_h, oi_m, oi_l;
    tw_add(w_ih, w_im, w_il, ti_h, ti_m, ti_l, oi_h, oi_m, oi_l);

    TW3_STORE6(v, base, or_h, oi_h, or_m, oi_m, or_l, oi_l);
  }
}

void aypx_tw3(real_t a, real_t *y, int nv1,
              real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);
  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  aypx_kernel_tw3<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);
  cudaDeviceSynchronize();
}

void aypx_tw3_pair(real_t a_h, real_t a_m, real_t a_l,
                   real_t *y, int nv1,
                   real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);
  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  aypx_kernel_tw3_pair<<<nbl, nth>>>(a_h, a_m, a_l, &y_dev[nv1], &x_dev[nv2], nin);
  cudaDeviceSynchronize();
}

//====================================================================
// QTW norm2: 3-word thread accumulator over QTW field.
// Per complex, compute (r^2 + i^2) where r = rh+rm+rl, i = ih+im+il via
// tw_mul on the TW reals, then tw_add into the accumulator.
__global__ void norm2_reduce_fused_kernel_tw3(real_t *red,
                                              const real_t *__restrict__ v1,
                                              int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sh = 0, sm = 0, sl = 0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      const int nin6 = nin / 6;
      for (int in6 = 0; in6 < nin6; ++in6)
      {
        int base = 6 * IDX2(nin6, in6, ist2);
        TW3_LOAD6(v, v1, base);

        // sq_r = (v_rh, v_rm, v_rl)^2
        real_t sqr_h, sqr_m, sqr_l;
        tw_mul(v_rh, v_rm, v_rl, v_rh, v_rm, v_rl, sqr_h, sqr_m, sqr_l);
        // sq_i = (v_ih, v_im, v_il)^2
        real_t sqi_h, sqi_m, sqi_l;
        tw_mul(v_ih, v_im, v_il, v_ih, v_im, v_il, sqi_h, sqi_m, sqi_l);
        // contribution = sq_r + sq_i
        real_t c_h, c_m, c_l;
        tw_add(sqr_h, sqr_m, sqr_l, sqi_h, sqi_m, sqi_l, c_h, c_m, c_l);
        // sum += contribution
        real_t ns_h, ns_m, ns_l;
        tw_add(sh, sm, sl, c_h, c_m, c_l, ns_h, ns_m, ns_l);
        sh = ns_h; sm = ns_m; sl = ns_l;
      }
    }
  }

  real_t *sh_h = sdata;
  real_t *sh_m = sdata +     blockDim.x;
  real_t *sh_l = sdata + 2 * blockDim.x;
  sh_h[tid] = sh;
  sh_m[tid] = sm;
  sh_l[tid] = sl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t bh = sh_h[tid+s], bm = sh_m[tid+s], bl = sh_l[tid+s];
      real_t ah = sh_h[tid],   am = sh_m[tid],   al = sh_l[tid];
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      sh_h[tid] = rh; sh_m[tid] = rm; sh_l[tid] = rl;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_h[tid], am = sh_m[tid], al = sh_l[tid];
#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      real_t bh = __shfl_down_sync(0xffffffff, ah, offset);
      real_t bm = __shfl_down_sync(0xffffffff, am, offset);
      real_t bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      ah = rh; am = rm; al = rl;
    }
    if (tid == 0)
    {
      red[blockIdx.x]                            = ah;
      red[blockIdx.x +     MAX_THREAD_PER_BLOCK] = am;
      red[blockIdx.x + 2 * MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
// dot_tw3: real-valued inner product, 3-word thread accumulator.
__global__ void dot_reduce_fused_kernel_tw3(real_t *red,
                                            const real_t *__restrict__ v1,
                                            const real_t *__restrict__ v2,
                                            int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sh = 0, sm = 0, sl = 0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      const int nin6 = nin / 6;
      for (int in6 = 0; in6 < nin6; ++in6)
      {
        int base = 6 * IDX2(nin6, in6, ist2);
        TW3_LOAD6(a, v1, base);
        TW3_LOAD6(b, v2, base);

        // v1.r * v2.r
        real_t pr_h, pr_m, pr_l;
        tw_mul(a_rh, a_rm, a_rl, b_rh, b_rm, b_rl, pr_h, pr_m, pr_l);
        // v1.i * v2.i
        real_t pi_h, pi_m, pi_l;
        tw_mul(a_ih, a_im, a_il, b_ih, b_im, b_il, pi_h, pi_m, pi_l);
        // c = pr + pi
        real_t c_h, c_m, c_l;
        tw_add(pr_h, pr_m, pr_l, pi_h, pi_m, pi_l, c_h, c_m, c_l);
        // sum += c
        real_t ns_h, ns_m, ns_l;
        tw_add(sh, sm, sl, c_h, c_m, c_l, ns_h, ns_m, ns_l);
        sh = ns_h; sm = ns_m; sl = ns_l;
      }
    }
  }

  real_t *sh_h = sdata;
  real_t *sh_m = sdata +     blockDim.x;
  real_t *sh_l = sdata + 2 * blockDim.x;
  sh_h[tid] = sh;
  sh_m[tid] = sm;
  sh_l[tid] = sl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t bh = sh_h[tid+s], bm = sh_m[tid+s], bl = sh_l[tid+s];
      real_t ah = sh_h[tid],   am = sh_m[tid],   al = sh_l[tid];
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      sh_h[tid] = rh; sh_m[tid] = rm; sh_l[tid] = rl;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah = sh_h[tid], am = sh_m[tid], al = sh_l[tid];
#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      real_t bh = __shfl_down_sync(0xffffffff, ah, offset);
      real_t bm = __shfl_down_sync(0xffffffff, am, offset);
      real_t bl = __shfl_down_sync(0xffffffff, al, offset);
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      ah = rh; am = rm; al = rl;
    }
    if (tid == 0)
    {
      red[blockIdx.x]                            = ah;
      red[blockIdx.x +     MAX_THREAD_PER_BLOCK] = am;
      red[blockIdx.x + 2 * MAX_THREAD_PER_BLOCK] = al;
    }
  }
}

//====================================================================
// dotc_tw3: complex inner product, 3-word accumulator for re and im.
__global__ void dotc_reduce_fused_kernel_tw3(real_t *red1, real_t *red2,
                                             const real_t *__restrict__ v1,
                                             const real_t *__restrict__ v2,
                                             int nin, int nvol, int nex)
{
  extern __shared__ real_t sdata[];
  const int tid = threadIdx.x;
  const int ist = blockIdx.x * blockDim.x + threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t arh = 0, arm = 0, arl = 0;
  real_t aih = 0, aim = 0, ail = 0;

  for (int idx = ist; idx < nvol; idx += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = idx + nvol * ex;
      const int nin6 = nin / 6;
      for (int in6 = 0; in6 < nin6; ++in6)
      {
        int base = 6 * IDX2(nin6, in6, ist2);
        TW3_LOAD6(a, v1, base);
        TW3_LOAD6(b, v2, base);

        // Real: v1.r*v2.r + v1.i*v2.i
        real_t pr_h, pr_m, pr_l;
        tw_mul(a_rh, a_rm, a_rl, b_rh, b_rm, b_rl, pr_h, pr_m, pr_l);
        real_t pi_h, pi_m, pi_l;
        tw_mul(a_ih, a_im, a_il, b_ih, b_im, b_il, pi_h, pi_m, pi_l);
        real_t cr_h, cr_m, cr_l;
        tw_add(pr_h, pr_m, pr_l, pi_h, pi_m, pi_l, cr_h, cr_m, cr_l);
        real_t nrh, nrm, nrl;
        tw_add(arh, arm, arl, cr_h, cr_m, cr_l, nrh, nrm, nrl);
        arh = nrh; arm = nrm; arl = nrl;

        // Imag: v1.r*v2.i - v1.i*v2.r
        real_t qr_h, qr_m, qr_l;
        tw_mul(a_rh, a_rm, a_rl, b_ih, b_im, b_il, qr_h, qr_m, qr_l);
        real_t qi_h, qi_m, qi_l;
        tw_mul(a_ih, a_im, a_il, b_rh, b_rm, b_rl, qi_h, qi_m, qi_l);
        real_t ci_h, ci_m, ci_l;
        tw_add(qr_h, qr_m, qr_l, -qi_h, -qi_m, -qi_l, ci_h, ci_m, ci_l);
        real_t nih, nim, nil;
        tw_add(aih, aim, ail, ci_h, ci_m, ci_l, nih, nim, nil);
        aih = nih; aim = nim; ail = nil;
      }
    }
  }

  real_t *shr_h = sdata;
  real_t *shr_m = sdata +     blockDim.x;
  real_t *shr_l = sdata + 2 * blockDim.x;
  real_t *shi_h = sdata + 3 * blockDim.x;
  real_t *shi_m = sdata + 4 * blockDim.x;
  real_t *shi_l = sdata + 5 * blockDim.x;
  shr_h[tid] = arh; shr_m[tid] = arm; shr_l[tid] = arl;
  shi_h[tid] = aih; shi_m[tid] = aim; shi_l[tid] = ail;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s >= WARP_LENGTH; s >>= 1)
  {
    if (tid < s)
    {
      real_t bh = shr_h[tid+s], bm = shr_m[tid+s], bl = shr_l[tid+s];
      real_t ah = shr_h[tid],   am = shr_m[tid],   al = shr_l[tid];
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      shr_h[tid] = rh; shr_m[tid] = rm; shr_l[tid] = rl;

      bh = shi_h[tid+s]; bm = shi_m[tid+s]; bl = shi_l[tid+s];
      ah = shi_h[tid];   am = shi_m[tid];   al = shi_l[tid];
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      shi_h[tid] = rh; shi_m[tid] = rm; shi_l[tid] = rl;
    }
    __syncthreads();
  }

  if (tid < WARP_LENGTH)
  {
    real_t ah_r = shr_h[tid], am_r = shr_m[tid], al_r = shr_l[tid];
    real_t ah_i = shi_h[tid], am_i = shi_m[tid], al_i = shi_l[tid];
#pragma unroll
    for (int offset = 16; offset > 0; offset >>= 1)
    {
      real_t bh, bm, bl, rh, rm, rl;
      bh = __shfl_down_sync(0xffffffff, ah_r, offset);
      bm = __shfl_down_sync(0xffffffff, am_r, offset);
      bl = __shfl_down_sync(0xffffffff, al_r, offset);
      tw_add(ah_r, am_r, al_r, bh, bm, bl, rh, rm, rl);
      ah_r = rh; am_r = rm; al_r = rl;
      bh = __shfl_down_sync(0xffffffff, ah_i, offset);
      bm = __shfl_down_sync(0xffffffff, am_i, offset);
      bl = __shfl_down_sync(0xffffffff, al_i, offset);
      tw_add(ah_i, am_i, al_i, bh, bm, bl, rh, rm, rl);
      ah_i = rh; am_i = rm; al_i = rl;
    }
    if (tid == 0)
    {
      red1[blockIdx.x]                            = ah_r;
      red1[blockIdx.x +     MAX_THREAD_PER_BLOCK] = am_r;
      red1[blockIdx.x + 2 * MAX_THREAD_PER_BLOCK] = al_r;
      red2[blockIdx.x]                            = ah_i;
      red2[blockIdx.x +     MAX_THREAD_PER_BLOCK] = am_i;
      red2[blockIdx.x + 2 * MAX_THREAD_PER_BLOCK] = al_i;
    }
  }
}

//====================================================================
// QTW final-pass reducer: aggregates per-block (h, m, l) into a single
// triple stored at red[0], red[1], red[2].
__global__ void reduce_kernel_multiblocks_tw3(real_t *red, int n)
{
  extern __shared__ real_t sdata[];
  real_t *sh_h = sdata;
  real_t *sh_m = sdata +     blockDim.x;
  real_t *sh_l = sdata + 2 * blockDim.x;

  const int tid = threadIdx.x;
  const int gridSize = blockDim.x * gridDim.x;

  real_t sh = 0, sm = 0, sl = 0;
  for (int i = tid; i < n; i += gridSize)
  {
    real_t ah = red[i];
    real_t am = red[i +     MAX_THREAD_PER_BLOCK];
    real_t al = red[i + 2 * MAX_THREAD_PER_BLOCK];
    real_t rh, rm, rl;
    tw_add(sh, sm, sl, ah, am, al, rh, rm, rl);
    sh = rh; sm = rm; sl = rl;
  }

  sh_h[tid] = sh; sh_m[tid] = sm; sh_l[tid] = sl;
  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
  {
    if (tid < s)
    {
      real_t bh = sh_h[tid+s], bm = sh_m[tid+s], bl = sh_l[tid+s];
      real_t ah = sh_h[tid],   am = sh_m[tid],   al = sh_l[tid];
      real_t rh, rm, rl;
      tw_add(ah, am, al, bh, bm, bl, rh, rm, rl);
      sh_h[tid] = rh; sh_m[tid] = rm; sh_l[tid] = rl;
    }
    __syncthreads();
  }

  if (tid == 0)
  {
    red[0] = sh_h[0];
    red[1] = sh_m[0];
    red[2] = sh_l[0];
  }
}

//====================================================================
// Host wrappers for QTW BLAS reductions.
void norm2_tw3(real_t *h, real_t *m, real_t *l,
               real_t *x1, real_t *red1, int nin, int nvol, int nex)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  norm2_reduce_fused_kernel_tw3<<<blocksPerGrid, threadsPerBlock, 3*sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
  reduce_kernel_multiblocks_tw3<<<1, threadsPerBlock, 3*sharedMemSize>>>(red1_dev, blocksPerGrid);

  cudaDeviceSynchronize();

  real_t triple[3] = {0, 0, 0};
  CHECK(cudaMemcpy(triple, &red1_dev[0], 3 * sizeof(real_t), cudaMemcpyDeviceToHost));
  *h = triple[0]; *m = triple[1]; *l = triple[2];
}

void dot_tw3(real_t *h, real_t *m, real_t *l,
             real_t *x1, real_t *x2, real_t *red1,
             int nin, int nvol, int nex)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *x2_dev   = (real_t *)dev_ptr(x2);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  dot_reduce_fused_kernel_tw3<<<blocksPerGrid, threadsPerBlock, 3*sharedMemSize>>>(red1_dev, x1_dev, x2_dev, nin, nvol, nex);
  reduce_kernel_multiblocks_tw3<<<1, threadsPerBlock, 3*sharedMemSize>>>(red1_dev, blocksPerGrid);

  cudaDeviceSynchronize();

  real_t triple[3] = {0, 0, 0};
  CHECK(cudaMemcpy(triple, &red1_dev[0], 3 * sizeof(real_t), cudaMemcpyDeviceToHost));
  *h = triple[0]; *m = triple[1]; *l = triple[2];
}

void dotc_tw3(real_t *ar_h, real_t *ar_m, real_t *ar_l,
              real_t *ai_h, real_t *ai_m, real_t *ai_l,
              real_t *x1, real_t *x2, real_t *red1, real_t *red2,
              int nin, int nvol, int nex)
{
  real_t *x1_dev   = (real_t *)dev_ptr(x1);
  real_t *x2_dev   = (real_t *)dev_ptr(x2);
  real_t *red1_dev = (real_t *)dev_ptr(red1);
  real_t *red2_dev = (real_t *)dev_ptr(red2);

  int threadsPerBlock = VECTOR_LENGTH;
  int blocksPerGrid   = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  dotc_reduce_fused_kernel_tw3<<<blocksPerGrid, threadsPerBlock, 6*sharedMemSize>>>(red1_dev, red2_dev, x1_dev, x2_dev, nin, nvol, nex);
  reduce_kernel_multiblocks_tw3<<<1, threadsPerBlock, 3*sharedMemSize>>>(red1_dev, blocksPerGrid);
  reduce_kernel_multiblocks_tw3<<<1, threadsPerBlock, 3*sharedMemSize>>>(red2_dev, blocksPerGrid);

  cudaDeviceSynchronize();

  real_t r1[3] = {0,0,0}, r2[3] = {0,0,0};
  CHECK(cudaMemcpy(r1, &red1_dev[0], 3 * sizeof(real_t), cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(r2, &red2_dev[0], 3 * sizeof(real_t), cudaMemcpyDeviceToHost));
  *ar_h = r1[0]; *ar_m = r1[1]; *ar_l = r1[2];
  *ai_h = r2[0]; *ai_m = r2[1]; *ai_l = r2[2];
}

//====================================================================
// QTW renormalize: re-establish |m| <= ulp(h)/2 and |l| <= ulp(m)/2 after
// an arithmetic chain that may have left the triple loosely renormalized.
__global__ void normalize_kernel_tw3(real_t *spinor, int nin, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site >= nvol) return;
  const int nin6 = nin / 6;
  for (int in6 = 0; in6 < nin6; ++in6)
  {
    int base = 6 * IDX2(nin6, in6, site);
    real_t rh = spinor[base + 0], ih = spinor[base + 1];
    real_t rm = spinor[base + 2], im = spinor[base + 3];
    real_t rl = spinor[base + 4], il = spinor[base + 5];
    Renormalize3(rh, rm, rl);
    Renormalize3(ih, im, il);
    spinor[base + 0] = rh; spinor[base + 1] = ih;
    spinor[base + 2] = rm; spinor[base + 3] = im;
    spinor[base + 4] = rl; spinor[base + 5] = il;
  }
}

void normalize_tw3(real_t *v, int nin, int nvol)
{
  real_t *v_dev = (real_t *)dev_ptr(v);
  int nth = VECTOR_LENGTH;
  int nbl = (nvol + nth - 1) / nth;
  normalize_kernel_tw3<<<nbl, nth>>>(v_dev, nin, nvol);
  cudaDeviceSynchronize();
}

#undef TW3_LOAD6
#undef TW3_STORE6

//====================================================================
#ifdef __CUDACC__
namespace {
__global__ void normalize_kernel_qdw_actual(real_t *spinor, int nin, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if (site >= nvol) return;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    real_t h_re = spinor[base + 0];
    real_t h_im = spinor[base + 1];
    real_t l_re = spinor[base + 2];
    real_t l_im = spinor[base + 3];

    real_t nh_re, nl_re, nh_im, nl_im;
    TwoSum_r(h_re, l_re, nh_re, nl_re);
    TwoSum_r(h_im, l_im, nh_im, nl_im);

    spinor[base + 0] = nh_re;
    spinor[base + 1] = nh_im;
    spinor[base + 2] = nl_re;
    spinor[base + 3] = nl_im;
  }
}
}
#endif

void normalize(real_t* v, int nin, int nvol, int mode)
{
  // mode == MWMode::DW (1) → QDW two-word renormalize.
  // mode == MWMode::TW (2) → QTW three-word renormalize.
  // anything else (incl. FP=0) → no-op.
#ifdef __CUDACC__
  if (mode == 1) {
    real_t *v_dev = (real_t *)dev_ptr(v);
    int nth = VECTOR_LENGTH;
    int nbl = (nvol + nth - 1) / nth;
    normalize_kernel_qdw_actual<<<nbl, nth>>>(v_dev, nin, nvol);
    cudaDeviceSynchronize();
  } else if (mode == 2) {
    normalize_tw3(v, nin, nvol);
  }
#endif
}

//============================================================END=====
