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
__global__ void copy_to_qdw_kernel(double4 *qdw, const real_t *std_v, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if(site >= nvol) return;

  for(int id=0; id<ND; ++id) {
    for(int ic=0; ic<NC; ++ic) {
       int qdw_idx = IDX2_QDW(ic, id, site);
       int r_idx = IDX2_SP_R(ic, id, site);
       int i_idx = IDX2_SP_I(ic, id, site);
       
       double4 v;
       v.x = std_v[r_idx];
       v.y = std_v[i_idx];
       v.z = 0.0;
       v.w = 0.0;
       qdw[qdw_idx] = v;
    }
  }
}

//====================================================================
__global__ void copy_from_qdw_kernel(real_t *std_v, const double4 *qdw, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if(site >= nvol) return;

  for(int id=0; id<ND; ++id) {
    for(int ic=0; ic<NC; ++ic) {
       int qdw_idx = IDX2_QDW(ic, id, site);
       int r_idx = IDX2_SP_R(ic, id, site);
       int i_idx = IDX2_SP_I(ic, id, site);
       
       double4 v = qdw[qdw_idx];
       std_v[r_idx] = v.x + v.z;
       std_v[i_idx] = v.y + v.w;
    }
  }
}

//====================================================================
void copy_to_qdw(real_t* qdw_v, const real_t* std_v, int nvol)
{
  double4* qdw_dev = (double4*)dev_ptr(qdw_v);
  real_t* std_dev  = (real_t*)dev_ptr(const_cast<real_t*>(std_v));

  int nth = VECTOR_LENGTH;
  int nbl = CEIL_NWP(nvol) / nth; 
  if (nbl == 0) nbl = 1;

  copy_to_qdw_kernel<<<nbl, nth>>>(qdw_dev, std_dev, nvol);
  cudaDeviceSynchronize();
}

//====================================================================
void copy_from_qdw(real_t* std_v, const real_t* qdw_v, int nvol)
{
  real_t* std_dev  = (real_t*)dev_ptr(std_v);
  double4* qdw_dev = (double4*)dev_ptr(const_cast<real_t*>(qdw_v));

  int nth = VECTOR_LENGTH;
  int nbl = CEIL_NWP(nvol) / nth; 
  if (nbl == 0) nbl = 1;

  copy_from_qdw_kernel<<<nbl, nth>>>(std_dev, qdw_dev, nvol);
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

  for (int in4 = 0; in4 < nin / 4; ++in4)
  {
    double wh_re = (double)w[IDX2(nin, 4 * in4 + 0, site)];
    double wh_im = (double)w[IDX2(nin, 4 * in4 + 1, site)];
    double wl_re = (double)w[IDX2(nin, 4 * in4 + 2, site)];
    double wl_im = (double)w[IDX2(nin, 4 * in4 + 3, site)];
    
    double vh_re = (double)v[IDX2(nin, 4 * in4 + 0, site)];
    double vh_im = (double)v[IDX2(nin, 4 * in4 + 1, site)];
    double vl_re = (double)v[IDX2(nin, 4 * in4 + 2, site)];
    double vl_im = (double)v[IDX2(nin, 4 * in4 + 3, site)];

    double th_re, tl_re;
    dw_scal(a_d, wh_re, wl_re, th_re, tl_re);
    double th_im, tl_im;
    dw_scal(a_d, wh_im, wl_im, th_im, tl_im);
    
    double out_h_re, out_l_re;
    dw_add(vh_re, vl_re, th_re, tl_re, out_h_re, out_l_re);
    double out_h_im, out_l_im;
    dw_add(vh_im, vl_im, th_im, tl_im, out_h_im, out_l_im);
    
    v[IDX2(nin, 4 * in4 + 0, site)] = (real_t)out_h_re;
    v[IDX2(nin, 4 * in4 + 1, site)] = (real_t)out_h_im;
    v[IDX2(nin, 4 * in4 + 2, site)] = (real_t)out_l_re;
    v[IDX2(nin, 4 * in4 + 3, site)] = (real_t)out_l_im;
  }
}

//====================================================================
void axpy(real_t *y, int nv1, real_t a,
          real_t *x, int nv2, int nin, int nvol, int use_qdw)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (use_qdw) {
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
          real_t *x, int nv2, int nin, int nvol, int use_qdw)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (use_qdw) {
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

    double th_re, tl_re;
    dw_scal(a_d, vh_re, vl_re, th_re, tl_re);
    double th_im, tl_im;
    dw_scal(a_d, vh_im, vl_im, th_im, tl_im);
    
    double out_h_re, out_l_re;
    dw_add(wh_re, wl_re, th_re, tl_re, out_h_re, out_l_re);
    double out_h_im, out_l_im;
    dw_add(wh_im, wl_im, th_im, tl_im, out_h_im, out_l_im);
    
    v[IDX2(nin, 4 * in4 + 0, site)] = (real_t)out_h_re;
    v[IDX2(nin, 4 * in4 + 1, site)] = (real_t)out_h_im;
    v[IDX2(nin, 4 * in4 + 2, site)] = (real_t)out_l_re;
    v[IDX2(nin, 4 * in4 + 3, site)] = (real_t)out_l_im;
  }
}

//====================================================================
void aypx(real_t a, real_t *y, int nv1,
          real_t *x, int nv2, int nin, int nvol, int use_qdw)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (use_qdw) {
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
          real_t *x, int nv2, int nin, int nvol, int use_qdw)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  if (use_qdw) {
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
      
      // QDW variables are loaded 4 floats contiguous (1 double4)
      for (int in4 = 0; in4 < nin / 4; ++in4)
      {
        double h_re = (double)v1[IDX2(nin, 4 * in4 + 0, ist2)];
        double h_im = (double)v1[IDX2(nin, 4 * in4 + 1, ist2)];
        double l_re = (double)v1[IDX2(nin, 4 * in4 + 2, ist2)];
        double l_im = (double)v1[IDX2(nin, 4 * in4 + 3, ist2)];
        
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
real_t norm2(real_t *x1, real_t *red1, int nin, int nvol, int nex, int use_qdw)
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

  if (use_qdw) {
    norm2_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
  } else {
    norm2_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
  }
  reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);

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
      for (int in4 = 0; in4 < nin / 4; ++in4)
      {
        double v1h_re = (double)v1[IDX2(nin, 4 * in4 + 0, ist2)];
        double v1h_im = (double)v1[IDX2(nin, 4 * in4 + 1, ist2)];
        double v1l_re = (double)v1[IDX2(nin, 4 * in4 + 2, ist2)];
        double v1l_im = (double)v1[IDX2(nin, 4 * in4 + 3, ist2)];
        
        double v2h_re = (double)v2[IDX2(nin, 4 * in4 + 0, ist2)];
        double v2h_im = (double)v2[IDX2(nin, 4 * in4 + 1, ist2)];
        double v2l_re = (double)v2[IDX2(nin, 4 * in4 + 2, ist2)];
        double v2l_im = (double)v2[IDX2(nin, 4 * in4 + 3, ist2)];

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
real_t dot(real_t *x1, real_t *x2, real_t *red1,
           int nin, int nvol, int nex, int use_qdw)
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

  if (use_qdw) {
    dot_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev,
                                                                                   x1_dev, x2_dev,
                                                                                   nin, nvol, nex);
  } else {
    dot_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev,
                                                                              x1_dev, x2_dev,
                                                                              nin, nvol, nex);
  }

  reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
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
  real_t *sh_r = sdata;
  real_t *sh_i = sdata + blockDim.x;

  const int nin2 = nin / 2;

  real_t ar = 0.0, ai = 0.0;
  real_t comp_r = 0.0, comp_i = 0.0;  // Kahan compensation

  for (int ist = idx; ist < nvol; ist += gridSize)
  {
    for (int ex = 0; ex < nex; ++ex)
    {
      int ist2 = ist + nvol * ex;

      for (int in4 = 0; in4 < nin / 4; ++in4)
      {
        double v1h_re = (double)v1[IDX2(nin, 4 * in4 + 0, ist2)];
        double v1h_im = (double)v1[IDX2(nin, 4 * in4 + 1, ist2)];
        double v1l_re = (double)v1[IDX2(nin, 4 * in4 + 2, ist2)];
        double v1l_im = (double)v1[IDX2(nin, 4 * in4 + 3, ist2)];
        
        double v2h_re = (double)v2[IDX2(nin, 4 * in4 + 0, ist2)];
        double v2h_im = (double)v2[IDX2(nin, 4 * in4 + 1, ist2)];
        double v2l_re = (double)v2[IDX2(nin, 4 * in4 + 2, ist2)];
        double v2l_im = (double)v2[IDX2(nin, 4 * in4 + 3, ist2)];

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
void dotc(real_t *ar, real_t *ai, real_t *x1, real_t *x2,
          real_t *red1, real_t *red2, int nin, int nvol, int nex, int use_qdw)
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

  if (use_qdw) {
    dotc_reduce_fused_kernel_qdw<<<blocksPerGrid, threadsPerBlock, 2 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                        x1_dev, x2_dev, nin, nvol, nex);
  } else {
    dotc_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, 2 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                    x1_dev, x2_dev, nin, nvol, nex);
  }

  reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red2_dev, blocksPerGrid);
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
__global__ void normalize_kernel_qdw_actual(double4 *spinor, int nvol)
{
  const int site = blockIdx.x * blockDim.x + threadIdx.x;
  if(site >= nvol) return;

  for (int id = 0; id < ND; ++id) {
    for (int ic = 0; ic < NC; ++ic) {
      int idx = IDX2_QDW(ic, id, site);
      double4 s = spinor[idx];
      double sr, er, si, ei;
      TwoSum(s.x, s.z, sr, er);
      TwoSum(s.y, s.w, si, ei);
      s.x = sr; s.z = er;
      s.y = si; s.w = ei;
      spinor[idx] = s;
    }
  }
}
}
#endif

void normalize(real_t* v, int nin, int nvol, int use_qdw)
{
  if (!use_qdw) return;
#ifdef __CUDACC__
  if (sizeof(real_t) == 8) {
    double4 *v_dev = (double4 *)dev_ptr(v);
    int nth = VECTOR_LENGTH;
    int nbl = (nvol + nth - 1) / nth;
    normalize_kernel_qdw_actual<<<nbl, nth>>>(v_dev, nvol);
    cudaDeviceSynchronize();
  }
#endif
}

//============================================================END=====
