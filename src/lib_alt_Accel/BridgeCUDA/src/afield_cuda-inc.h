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
  double a_d = (double)a;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    double wh_re = (double)w[base + 0];
    double wh_im = (double)w[base + 1];
    double wl_re = (double)w[base + 2];
    double wl_im = (double)w[base + 3];

    double vh_re = (double)v[base + 0];
    double vh_im = (double)v[base + 1];
    double vl_re = (double)v[base + 2];
    double vl_im = (double)v[base + 3];

    double th_re, tl_re;
    dw_scal(a_d, wh_re, wl_re, th_re, tl_re);
    double th_im, tl_im;
    dw_scal(a_d, wh_im, wl_im, th_im, tl_im);

    double out_h_re, out_l_re;
    dw_add(vh_re, vl_re, th_re, tl_re, out_h_re, out_l_re);
    double out_h_im, out_l_im;
    dw_add(vh_im, vl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = (real_t)out_h_re;
    v[base + 1] = (real_t)out_h_im;
    v[base + 2] = (real_t)out_l_re;
    v[base + 3] = (real_t)out_l_im;
  }
}

//====================================================================
void axpy(real_t *y, int nv1, real_t a,
          real_t *x, int nv2, int nin, int nvol, int mode)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (mode == 1) {
    axpy_kernel_qdw<<<nbl, nth>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
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
  double a_d = (double)a;

  const int nin4 = nin / 4;
  for (int in4 = 0; in4 < nin4; ++in4)
  {
    int base = 4 * IDX2(nin4, in4, site);
    double vh_re = (double)v[base + 0];
    double vh_im = (double)v[base + 1];
    double vl_re = (double)v[base + 2];
    double vl_im = (double)v[base + 3];

    double wh_re = (double)w[base + 0];
    double wh_im = (double)w[base + 1];
    double wl_re = (double)w[base + 2];
    double wl_im = (double)w[base + 3];

    double th_re, tl_re;
    dw_scal(a_d, vh_re, vl_re, th_re, tl_re);
    double th_im, tl_im;
    dw_scal(a_d, vh_im, vl_im, th_im, tl_im);

    double out_h_re, out_l_re;
    dw_add(wh_re, wl_re, th_re, tl_re, out_h_re, out_l_re);
    double out_h_im, out_l_im;
    dw_add(wh_im, wl_im, th_im, tl_im, out_h_im, out_l_im);

    v[base + 0] = (real_t)out_h_re;
    v[base + 1] = (real_t)out_h_im;
    v[base + 2] = (real_t)out_l_re;
    v[base + 3] = (real_t)out_l_im;
  }
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
  double ar_d = (double)ar;
  double ai_d = (double)ai;

  for (int in4 = 0; in4 < nin / 4; ++in4)
  {
    double vh_re = (double)v[IDX2(nin, 4 * in4 + 0, site)];
    double vh_im = (double)v[IDX2(nin, 4 * in4 + 1, site)];
    double vl_re = (double)v[IDX2(nin, 4 * in4 + 2, site)];
    double vl_im = (double)v[IDX2(nin, 4 * in4 + 3, site)];
    
    double wh_re = (double)w[IDX2(nin, 4 * in4 + 0, site)];
    double wh_im = (double)w[IDX2(nin, 4 * in4 + 1, site)];
    double wl_re = (double)w[IDX2(nin, 4 * in4 + 2, site)];
    double wl_im = (double)w[IDX2(nin, 4 * in4 + 3, site)];

    // real part: ar*vr - ai*vi
    double p1_r, e1_r; dw_scal(ar_d, vh_re, vl_re, p1_r, e1_r);
    double p2_r, e2_r; dw_scal(ai_d, vh_im, vl_im, p2_r, e2_r);
    double r_h, r_l; dw_add(p1_r, e1_r, -p2_r, -e2_r, r_h, r_l);
    
    // imag part: ai*vr + ar*vi
    double p1_i, e1_i; dw_scal(ai_d, vh_re, vl_re, p1_i, e1_i);
    double p2_i, e2_i; dw_scal(ar_d, vh_im, vl_im, p2_i, e2_i);
    double i_h, i_l; dw_add(p1_i, e1_i, p2_i, e2_i, i_h, i_l);
    
    // add w: v = a*v + w
    double out_h_re, out_l_re; dw_add(wh_re, wl_re, r_h, r_l, out_h_re, out_l_re);
    double out_h_im, out_l_im; dw_add(wh_im, wl_im, i_h, i_l, out_h_im, out_l_im);
    
    v[IDX2(nin, 4 * in4 + 0, site)] = (real_t)out_h_re;
    v[IDX2(nin, 4 * in4 + 1, site)] = (real_t)out_h_im;
    v[IDX2(nin, 4 * in4 + 2, site)] = (real_t)out_l_re;
    v[IDX2(nin, 4 * in4 + 3, site)] = (real_t)out_l_im;
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
        double h_re = (double)v1[base + 0];
        double h_im = (double)v1[base + 1];
        double l_re = (double)v1[base + 2];
        double l_im = (double)v1[base + 3];
        
        // Exact real square
        double p1, e1;
        TwoProd(h_re, h_re, p1, e1);
        double lo = 2.0 * h_re * l_re + e1;
        
        double dh, dl;
        double sum_d = (double)sum;
        TwoSum(sum_d, p1, dh, dl);
        double temp_sum = dh;
        double temp_comp = (double)comp + dl + lo;
        double sum_out, comp_out;
        TwoSum(temp_sum, temp_comp, sum_out, comp_out);
        sum = (real_t)sum_out;
        comp = (real_t)comp_out;
        
        // Exact imaginary square
        TwoProd(h_im, h_im, p1, e1);
        lo = 2.0 * h_im * l_im + e1;
        
        sum_d = (double)sum;
        TwoSum(sum_d, p1, dh, dl);
        temp_sum = dh;
        temp_comp = (double)comp + dl + lo;
        TwoSum(temp_sum, temp_comp, sum_out, comp_out);
        sum = (real_t)sum_out;
        comp = (real_t)comp_out;
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
        double v1h_re = (double)v1[base + 0];
        double v1h_im = (double)v1[base + 1];
        double v1l_re = (double)v1[base + 2];
        double v1l_im = (double)v1[base + 3];

        double v2h_re = (double)v2[base + 0];
        double v2h_im = (double)v2[base + 1];
        double v2l_re = (double)v2[base + 2];
        double v2l_im = (double)v2[base + 3];

        // Real dot product part: v1.re * v2.re + v1.im * v2.im
        double p1_r, e1_r, p2_r, e2_r;
        TwoProd(v1h_re, v2h_re, p1_r, e1_r);
        TwoProd(v1h_im, v2h_im, p2_r, e2_r);
        
        double dh_r, dl_r;
        TwoSum(p1_r, p2_r, dh_r, dl_r);
        
        double lo_r = v1h_re * v2l_re + v1l_re * v2h_re + v1h_im * v2l_im + v1l_im * v2h_im + e1_r + e2_r + dl_r;
        
        double dh, dl;
        double at_d = (double)at;
        TwoSum(at_d, dh_r, dh, dl);
        double temp_sum = dh;
        double temp_comp = (double)comp + dl + lo_r;
        double at_out, comp_out;
        TwoSum(temp_sum, temp_comp, at_out, comp_out);
        at = (real_t)at_out;
        comp = (real_t)comp_out;
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
        double v1h_re = (double)v1[base + 0];
        double v1h_im = (double)v1[base + 1];
        double v1l_re = (double)v1[base + 2];
        double v1l_im = (double)v1[base + 3];

        double v2h_re = (double)v2[base + 0];
        double v2h_im = (double)v2[base + 1];
        double v2l_re = (double)v2[base + 2];
        double v2l_im = (double)v2[base + 3];

        // Real part: v1.re * v2.re + v1.im * v2.im
        double p1_r, e1_r, p2_r, e2_r;
        TwoProd(v1h_re, v2h_re, p1_r, e1_r);
        TwoProd(v1h_im, v2h_im, p2_r, e2_r);
        
        double dh_r, dl_r;
        TwoSum(p1_r, p2_r, dh_r, dl_r);
        
        double lo_r = v1h_re * v2l_re + v1l_re * v2h_re + v1h_im * v2l_im + v1l_im * v2h_im + e1_r + e2_r + dl_r;
        
        double dh, dl;
        double ar_d = (double)ar;
        TwoSum(ar_d, dh_r, dh, dl);
        double temp_sum_r = dh;
        double temp_comp_r = (double)comp_r + dl + lo_r;
        double ar_out, comp_r_out;
        TwoSum(temp_sum_r, temp_comp_r, ar_out, comp_r_out);
        ar = (real_t)ar_out;
        comp_r = (real_t)comp_r_out;

        // Imag part: v1.re * v2.im - v1.im * v2.re
        double p1_i, e1_i, p2_i, e2_i;
        TwoProd(v1h_re, v2h_im, p1_i, e1_i);
        TwoProd(-v1h_im, v2h_re, p2_i, e2_i);
        
        double dh_i, dl_i;
        TwoSum(p1_i, p2_i, dh_i, dl_i);
        
        double lo_i = v1h_re * v2l_im + v1l_re * v2h_im - v1h_im * v2l_re - v1l_im * v2h_re + e1_i + e2_i + dl_i;
        
        double ai_d = (double)ai;
        TwoSum(ai_d, dh_i, dh, dl);
        double temp_sum_i = dh;
        double temp_comp_i = (double)comp_i + dl + lo_i;
        double ai_out, comp_i_out;
        TwoSum(temp_sum_i, temp_comp_i, ai_out, comp_i_out);
        ai = (real_t)ai_out;
        comp_i = (real_t)comp_i_out;
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
  if (mode != 1) return;
#ifdef __CUDACC__
  real_t *v_dev = (real_t *)dev_ptr(v);
  int nth = VECTOR_LENGTH;
  int nbl = (nvol + nth - 1) / nth;
  normalize_kernel_qdw_actual<<<nbl, nth>>>(v_dev, nin, nvol);
  cudaDeviceSynchronize();
#endif
}

//============================================================END=====
