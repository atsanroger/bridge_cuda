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
void axpy(real_t *y, int nv1, real_t a,
          real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  size_t sharedMemSize = nth * NWP * sizeof(real_t);

  // axpy_kernel<<<nbl,nth>>>(&y_dev[nv1], a, &x_dev[nv2], nin);
  axpy_kernel_opt<<<nbl, nth, sharedMemSize>>>(&y_dev[nv1], a, &x_dev[nv2], nin);

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
void axpy(real_t *y, int nv1, real_t ar, real_t ai,
          real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  axpy_kernel<<<nbl, nth>>>(&y_dev[nv1], ar, ai, &x_dev[nv2], nin);

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
void aypx(real_t a, real_t *y, int nv1,
          real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  aypx_kernel<<<nbl, nth>>>(a, &y_dev[nv1], &x_dev[nv2], nin);

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
void aypx(real_t ar, real_t ai, real_t *y, int nv1,
          real_t *x, int nv2, int nin, int nvol)
{
  real_t *y_dev = (real_t *)dev_ptr(y);
  real_t *x_dev = (real_t *)dev_ptr(x);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  aypx_kernel<<<nbl, nth>>>(ar, ai, &y_dev[nv1], &x_dev[nv2], nin);

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
real_t norm2(real_t *x1, real_t *red1, int nin, int nvol, int nex)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  // norm2_kernel<<<nbl,nth>>>(red1_dev, x1_dev, nin, nvol, nex);
  // cudaDeviceSynchronize();

  int nth2 = VECTOR_LENGTH;
  //  int nth2 = 1;
  int nbl2 = 1;
  // reduce_kernel<<<nbl2,nth2>>>(red1_dev, nvol);

  int threadsPerBlock = nth2;
  int blocksPerGrid = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  // reduce_kernel_multiblocks<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, nvol);
  // reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  norm2_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, x1_dev, nin, nvol, nex);
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
real_t dot(real_t *x1, real_t *x2, real_t *red1,
           int nin, int nvol, int nex)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *x2_dev = (real_t *)dev_ptr(x2);

  real_t *red1_dev = (real_t *)dev_ptr(red1);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;
  // dot_kernel<<<nbl,nth>>>(red1_dev, x1_dev, x2_dev, nin, nvol, nex);

  // cudaDeviceSynchronize();

  int nth2 = VECTOR_LENGTH;
  int nbl2 = 1;
  int nSM = getNumSMs();
  int threadsPerBlock = nth2;
  int blocksPerGrid = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  dot_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev,
                                                                             x1_dev, x2_dev,
                                                                             nin, nvol, nex);

  // reduce_kernel_multiblocks<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, nvol);
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
void dotc(real_t *ar, real_t *ai, real_t *x1, real_t *x2,
          real_t *red1, real_t *red2, int nin, int nvol, int nex)
{
  real_t *x1_dev = (real_t *)dev_ptr(x1);
  real_t *x2_dev = (real_t *)dev_ptr(x2);

  real_t *red1_dev = (real_t *)dev_ptr(red1);
  real_t *red2_dev = (real_t *)dev_ptr(red2);

  int nth = VECTOR_LENGTH;
  int nbl = nvol / nth;

  // dotc_kernel<<<nbl,nth>>>(red1_dev, red2_dev, x1_dev, x2_dev,
  //                          nin, nvol, nex);
  // cudaDeviceSynchronize();

  int nth2 = VECTOR_LENGTH;
  //  int nth2 = 1;
  int nbl2 = 1;
  // reduce_kernel<<<nbl2,nth2>>>(red1_dev, nvol);
  // reduce_kernel<<<nbl2,nth2>>>(red2_dev, nvol);
  int threadsPerBlock  = nth2;
  int blocksPerGrid    = min((nvol + threadsPerBlock - 1) / threadsPerBlock, MAX_THREAD_PER_BLOCK);
  size_t sharedMemSize = threadsPerBlock * sizeof(real_t);

  dotc_reduce_fused_kernel<<<blocksPerGrid, threadsPerBlock, 2 * sharedMemSize>>>(red1_dev, red2_dev,
                                                                                  x1_dev, x2_dev, nin, nvol, nex);
  // reduce_kernel_multiblocks<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red1_dev, nvol);
  reduce_kernel_multiblocks<<<nbl2, threadsPerBlock, sharedMemSize>>>(red1_dev, blocksPerGrid);
  // reduce_kernel_multiblocks<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(red2_dev, nvol);
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

//============================================================END=====
