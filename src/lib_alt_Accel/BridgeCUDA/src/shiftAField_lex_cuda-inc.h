/*!
      @file    shiftAField_lex_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2535 $
*/

//====================================================================
__global__
void shift_lex_xp1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = 0;
  int ist = ix + Nx * iyzt;

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    buf[in + nin * iyzt] = bc2 * v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_xp1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 0;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_xp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_xp2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = Nx-1;
  int ist = ix + Nx * iyzt;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = buf[in + nin * iyzt];
  }

}

//====================================================================
void shift_lex_xp2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 0;
  int bc1 = bc[idir];
  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_xp2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_xpb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix   = ist % Nx;
  int iyzt = ist / Nx;
  int nei  = ((ix+1) % Nx) + Nx * iyzt;

  real_t bc2 = 1.0;
  if(ix == Nx-1) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_xpb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 0;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_xpb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_xm1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = Nx-1;
  int ist = ix + Nx * iyzt;

  for(int in = 0; in < nin; ++in){
    buf[in + nin * iyzt] = v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_xm1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 0;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_xm1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_xm2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int iyzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = 0;
  int ist = ix + Nx * iyzt;

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * buf[in + nin * iyzt];
  }

}

//====================================================================
void shift_lex_xm2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 0;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_xm2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_xmb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix   = ist % Nx;
  int iyzt = ist/Nx;

  int nei  = ix-1 + Nx * iyzt;
  if(ix == 0) nei = Nx-1 + Nx * iyzt;

  real_t bc2 = 1.0;
  if(ix == 0) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_xmb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 0;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_xmb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_yp1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy   = 0;
  int ist  = ix + Nx * (iy + Ny * izt);

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
     buf[in + nin * ixzt] = bc2 * v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_yp1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 1;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_yp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_yp2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy   = Ny-1;
  int ist  = ix + Nx * (iy + Ny * izt);

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = buf[in + nin * ixzt];
  }

}

//====================================================================
void shift_lex_yp2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 1;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_yp2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_ypb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ist % Nx;
  int iy  = (ist/Nx) % Ny;
  int izt = ist/(Nx*Ny);
  int iy2 = (iy+1) % Ny;
  int nei = ix + Nx * (iy2 + Ny * izt);

  real_t bc2 = 1.0;
  if(iy == Ny-1) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_ypb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 1;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_ypb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_ym1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy  = Ny-1;
  int ist = ix + Nx * (iy + Ny*izt);

  for(int in = 0; in < nin; ++in){
    buf[in + nin * ixzt] = v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_ym1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 1;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_ym1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_ym2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixzt = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ixzt % Nx;
  int izt = ixzt / Nx;
  int iy = 0;
  int ist  = ix + Nx*(iy + Ny*izt);

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * buf[in + nin*ixzt];
  }

}

//====================================================================
void shift_lex_ym2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 1;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_ym2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_ymb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int ix  = ist % Nx;
  int iy  = (ist/Nx) % Ny;
  int izt = ist/(Nx*Ny);
  int iy2 = (iy-1+Ny) % Ny;
  int nei = ix + Nx * (iy2 + Ny * izt);

  real_t bc2 = 1.0;
  if(iy == 0) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_ymb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 1;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_ymb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_zp1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy = Nx * Ny;
  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz  = 0;
  int ist = ixy + Nxy * (iz + Nz * it);

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    buf[in + nin * ixyt] = bc2 * v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_zp1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 2;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_zp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_zp2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy = Nx * Ny;
  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz  = Nz-1;
  int ist = ixy + Nxy * (iz + Nz * it);

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = buf[in + nin * ixyt];
  }

}

//====================================================================
void shift_lex_zp2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 2;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_zp2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_zpb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy  = Nx * Ny;
  int ixy = ist % Nxy;
  int iz  = (ist/Nxy) % Nz;
  int it  = ist/(Nxy*Nz);
  int iz2 = (iz+1) % Nz;
  int nei = ixy + Nxy * (iz2 + Nz * it);

  real_t bc2 = 1.0;
  if(iz == Nz-1) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_zpb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 2;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_zpb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_zm1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy = Nx * Ny;
  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz  = Nz-1;
  int ist = ixy + Nxy * (iz + Nz * it);

  for(int in = 0; in < nin; ++in){
    buf[in + nin * ixyt] = v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_zm1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 2;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_zm1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_zm2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyt = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy = Nx * Ny;
  int ixy = ixyt % Nxy;
  int it  = ixyt / Nxy;
  int iz  = 0;
  int ist = ixy + Nxy * (iz + Nz * it);

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * buf[in + nin*ixyt];
  }

}

//====================================================================
void shift_lex_zm2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 2;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_zm2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//==================================================================== 
__global__
void shift_lex_zmb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxy = Nx * Ny;
  int ixy = ist % Nxy;
  int iz  = (ist/Nxy) % Nz;
  int it  = ist/(Nxy*Nz);
  int iz2 = (iz-1+Nz) % Nz;
  int nei = ixy + Nxy * (iz2 + Nz * it);

  real_t bc2 = 1.0;
  if(iz == 0) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_zmb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 2;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_zmb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_tp1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int it   = 0;
  int ist  = ixyz + Nxyz * it;

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    buf[in + nin * ixyz] = bc2 * v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_tp1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 3;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_tp1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_tp2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int it   = Nt-1;
  int ist  = ixyz + Nxyz * it;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = buf[in + nin * ixyz];
  }

}

//====================================================================
void shift_lex_tp2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 3;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_tp2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_tpb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int ixyz = ist % Nxyz;
  int it  = ist/Nxyz;
  int it2 = (it+1) % Nt;
  int nei = ixyz + Nxyz * it2;

  real_t bc2 = 1.0;
  if(it == Nt-1) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_tpb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 3;
  int bc1  = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_tpb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_tm1_kernel(real_t* buf, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int it   = Nt-1;
  int ist  = ixyz + Nxyz * it;

  for(int in = 0; in < nin; ++in){
    buf[in + nin * ixyz] = v1[IDX2(nin, in, ist)];
  }

}

//====================================================================
void shift_lex_tm1(real_t* buf, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* buf_dev = (real_t*)dev_ptr(buf);
  real_t* v1_dev  = (real_t*)dev_ptr(v1);

  int idir = 3;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_tm1_kernel<<<nbl,nth>>>(buf_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================
__global__
void shift_lex_tm2_kernel(real_t* v2, real_t* buf,  int nin, 
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ixyz = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int it   = 0;
  int ist  = ixyz + Nxyz * it;

  real_t bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * buf[in + nin * ixyz];
  }

}

//====================================================================
void shift_lex_tm2(real_t* v2, real_t* buf,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev  = (real_t*)dev_ptr(v2);
  real_t* buf_dev = (real_t*)dev_ptr(buf);
 
  int idir = 3;
  int bc1 = bc[idir];

  int Nbd = Nst/Nsize[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nbd/nth;
  shift_lex_tm2_kernel<<<nbl,nth>>>(v2_dev, buf_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//==================================================================== 
__global__
void shift_lex_tmb_kernel(real_t* v2, real_t* v1, int nin,
                          int Nx, int Ny, int Nz, int Nt, int bc1)
{
  int ist = blockIdx.x * blockDim.x + threadIdx.x;

  int Nxyz = Nx * Ny * Nz;
  int ixyz = ist % Nxyz;
  int it   = ist/Nxyz;
  int it2  = (it-1+Nt) % Nt;
  int nei  = ixyz + Nxyz * it2;

  real_t bc2 = 1.0;
  if(it == 0) bc2 = bc1;

  for(int in = 0; in < nin; ++in){
    v2[IDX2(nin, in, ist)] = bc2 * v1[IDX2(nin, in, nei)];
  }

}

//====================================================================
void shift_lex_tmb(real_t* v2, real_t* v1,
                   int nin, int *Nsize, int *bc)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  real_t* v2_dev = (real_t*)dev_ptr(v2);
  real_t* v1_dev = (real_t*)dev_ptr(v1);
 
  int idir = 3;
  int bc1 = bc[idir];

  int nth = VECTOR_LENGTH;
  int nbl = Nst/nth;
  shift_lex_tmb_kernel<<<nbl,nth>>>(v2_dev, v1_dev,
                                    nin, Nx, Ny, Nz, Nt, bc1);

}

//====================================================================

