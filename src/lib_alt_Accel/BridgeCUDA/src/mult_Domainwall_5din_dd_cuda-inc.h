/*!
      @file    mult_Doainwall_5din_dd_cuda-inc.h
      @brief
      @author  Wei-Lun Chen (wlchen)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2024-03-12 10:53:53 #$
      @version $LastChangedRevision: 2598 $
*/

#ifndef MULT_DOMAINWALL_5DIN_DD_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_DD_ACC_INCLUDED

//====================================================================
__global__ void mult_domainwall_5din_dd_5dir_dirac_dev(
                   real_t * vp, real_t * yp, real_t * wp,
                   real_t mq, real_t M0,
                   int Ns,
                   int Nx, int Ny, int Nz, int Nt,
                   int Bx, int By, int Bz, int Bt,
                   int ieo, real_t alpha)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  //int NBt  = Nt/Bt;
  //int Nblock = NBx * NBy * NBz * NBt;

  // Parallelize over (site, color-complex): one thread per (site2, ivc), ivc in
  // [0,NVC). The s5 loop stays inside the thread. The Wilson-5d term is color
  // diagonal, so each ivc is independent (it only couples the 4 Dirac spins at
  // the same ivc) -> NVC x more parallelism than per-site. (cf. eo 5dir kernel.)
  int gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Moebius coefficients from __constant__ memory (uploaded by the operator's
  // set_precond_parameters via setDomainwallConstantBC). Avoids the per-call
  // cudaMalloc/memcpy + shared-memory staging the old code used.
  const real_t* b = const_b;
  const real_t* c = const_c;

  if (gid < Nst * NVC) {

    int ivc   = gid % NVC;
    int site2 = gid / NVC;

    int bsite = site2 % Bsize;
    int block = site2 / Bsize;

    int ibx = block % NBx;
    int iby = (block/NBx) % NBy;
    int ibz = (block/(NBx * NBy)) % NBz;
    int ibt = block/(NBx * NBy * NBz);

    int ieo_skip = 1 - ieo;
    int jeo = (ibx + iby + ibz + ibt) % 2;
    if(jeo != ieo_skip){

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite / Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt / By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt / Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      for (int is = 0; is < Ns; ++is) {

        int is_up = (is+1) % Ns;
        real_t Fup = 0.5 * alpha;
        if (is == Ns-1) Fup = -0.5 * mq;

        real_t wu1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
        real_t wu2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
        real_t wu3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
        real_t wu4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];

        real_t vt1 = Fup * (wu1 - wu3);
        real_t vt2 = Fup * (wu2 - wu4);
        real_t vt3 = Fup * (wu3 - wu1);
        real_t vt4 = Fup * (wu4 - wu2);

        int is_dn = (is-1 + Ns) % Ns;
        real_t Fdn = 0.5 * alpha;
        if (is == 0) Fdn = -0.5 * mq;

        real_t wd1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
        real_t wd2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
        real_t wd3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
        real_t wd4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];

        vt1 += Fdn * (wd1 + wd3);
        vt2 += Fdn * (wd2 + wd4);
        vt3 += Fdn * (wd3 + wd1);
        vt4 += Fdn * (wd4 + wd2);

        real_t B_is = b[is] * (4.0 - M0) + 1.0;
        real_t C_is = c[is] * (4.0 - M0) - 1.0;

        real_t w1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        real_t w2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        real_t w3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        real_t w4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];

        // alpha congruence on the self (diagonal-in-s) block: bulk -> alpha*w,
        // s=0 / s=Ns-1 corners -> chiral (gamma5) mix 0.5(1+alpha) +/- 0.5(-+1+alpha).
        real_t st1, st2, st3, st4;
        if (is == 0) {
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (-1.0 + alpha);
          st1 = f1 * w1 + f2 * w3;  st2 = f1 * w2 + f2 * w4;
          st3 = f1 * w3 + f2 * w1;  st4 = f1 * w4 + f2 * w2;
        } else if (is == Ns - 1) {
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (1.0 - alpha);
          st1 = f1 * w1 + f2 * w3;  st2 = f1 * w2 + f2 * w4;
          st3 = f1 * w3 + f2 * w1;  st4 = f1 * w4 + f2 * w2;
        } else {
          st1 = alpha * w1;  st2 = alpha * w2;  st3 = alpha * w3;  st4 = alpha * w4;
        }

        vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = B_is * st1 + C_is * vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = B_is * st2 + C_is * vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = B_is * st3 + C_is * vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = B_is * st4 + C_is * vt4;

        yp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = -0.5 * (b[is]*st1 + c[is]*vt1);
        yp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = -0.5 * (b[is]*st2 + c[is]*vt2);
        yp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = -0.5 * (b[is]*st3 + c[is]*vt3);
        yp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = -0.5 * (b[is]*st4 + c[is]*vt4);
      }
    }
  }
}

//====================================================================
void mult_domainwall_5din_dd_5dir_dirac(
                   real_t *vp, real_t *yp, real_t *wp,
                   real_t mq, real_t M0,
                   int Ns,
                   real_t *b, real_t *c,
                   int *Nsize, int *block_size, int ieo, real_t alpha)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  int Nst  = Nx * Ny * Nz * Nt;

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  // b, c live in __constant__ memory (uploaded at set_config); kept as params
  // only for call-site ABI compatibility.
  (void)b; (void)c;

  // Define block and grid sizes for CUDA kernel: one thread per (site, color)
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst * NVC + blockSize - 1)/ blockSize;

  mult_domainwall_5din_dd_5dir_dirac_dev<<<gridSize, blockSize>>>(
                   vp_dev, yp_dev, wp_dev,
                   mq, M0, Ns,
                   Nx, Ny, Nz, Nt,
                   Bx, By, Bz, Bt,
                   ieo, alpha);

  //Synchronize to ensure completion of kernel execution
  dw5din_kernel_sync();
}

//====================================================================
__global__ void mult_domainwall_5din_dd_5dirdag_dirac_dev(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns,
       int Nx, int Ny, int Nz, int Nt,
       int Bx, int By, int Bz, int Bt,
       int ieo, real_t alpha)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;
  // size of block
  int Bsize = Bx * By * Bz * Bt;
  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  int NBt  = Nt/Bt;
  int Nblock = NBx * NBy * NBz * NBt;

  // Moebius coefficients from __constant__ memory (see dd_5dir_dirac_dev).
  const real_t* b = const_b;
  const real_t* c = const_c;

  // One thread per (site, color); see dd_5dir_dirac_dev.
  int gid = blockIdx.x * blockDim.x + threadIdx.x;

  if( gid < Nst * NVC ){

    int ivc   = gid % NVC;
    int site2 = gid / NVC;

    int bsite = site2 % Bsize;
    int block = site2/Bsize;
    int ibx = block % NBx;
    int iby = (block/NBx) % NBy;
    int ibz = (block/(NBx * NBy)) % NBz;
    int ibt = block/(NBx * NBy * NBz);

    int ieo_skip = 1 - ieo;
    int jeo = (ibx + iby + ibz + ibt) % 2;
    if(jeo != ieo_skip){

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite / Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt / By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt / Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      // color-parallel: one thread per (site, ivc); s5 loop inside. The 5d-dag
      // term is color diagonal so each ivc is independent (cf. dd_5dir_dirac_dev).
      for (int is = 0; is < Ns; ++is) {

        real_t B1 = b[is] * (4.0 - M0) + 1.0;
        real_t a1 = -0.5 * b[is];

        real_t wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        real_t wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        real_t wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        real_t wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];
        real_t yt1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
        real_t yt2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
        real_t yt3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
        real_t yt4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];

        // self block: S = B1*w + a1*gamma5(y); then alpha congruence (bulk
        // alpha*S; s=0/Ns-1 corners chirally mixed, same as dd_5dir_dirac_dev).
        real_t s1 = B1 * wt1 + a1 * yt3;
        real_t s2 = B1 * wt2 + a1 * yt4;
        real_t s3 = B1 * wt3 + a1 * yt1;
        real_t s4 = B1 * wt4 + a1 * yt2;
        real_t vt1, vt2, vt3, vt4;
        if (is == 0) {
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (-1.0 + alpha);
          vt1 = f1 * s1 + f2 * s3;  vt2 = f1 * s2 + f2 * s4;
          vt3 = f1 * s3 + f2 * s1;  vt4 = f1 * s4 + f2 * s2;
        } else if (is == Ns - 1) {
          real_t f1 = 0.5 * (1.0 + alpha), f2 = 0.5 * (1.0 - alpha);
          vt1 = f1 * s1 + f2 * s3;  vt2 = f1 * s2 + f2 * s4;
          vt3 = f1 * s3 + f2 * s1;  vt4 = f1 * s4 + f2 * s2;
        } else {
          vt1 = alpha * s1;  vt2 = alpha * s2;  vt3 = alpha * s3;  vt4 = alpha * s4;
        }

        int is_up = (is+1) % Ns;
        real_t C1 = c[is_up] * (4.0 - M0) - 1.0;
        real_t aup = -0.5 * c[is_up];

        real_t wu1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
        real_t wu2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
        real_t wu3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
        real_t wu4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];
        real_t yu1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
        real_t yu2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
        real_t yu3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
        real_t yu4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];

        real_t xu1 = C1 * wu1 + aup * yu3;
        real_t xu2 = C1 * wu2 + aup * yu4;
        real_t xu3 = C1 * wu3 + aup * yu1;
        real_t xu4 = C1 * wu4 + aup * yu2;

        real_t Fup = 0.5 * alpha;
        if (is == Ns-1) Fup = -0.5 * mq;

        vt1 += Fup * (xu1 + xu3);
        vt2 += Fup * (xu2 + xu4);
        vt3 += Fup * (xu3 + xu1);
        vt4 += Fup * (xu4 + xu2);

        int is_dn = (is-1 + Ns) % Ns;
        real_t C2 = c[is_dn] * (4.0 - M0) - 1.0;
        real_t adn = -0.5 * c[is_dn];

        real_t wd1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
        real_t wd2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
        real_t wd3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
        real_t wd4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];
        real_t yd1 = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
        real_t yd2 = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
        real_t yd3 = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
        real_t yd4 = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];

        real_t xd1 = C2 * wd1 + adn * yd3;
        real_t xd2 = C2 * wd2 + adn * yd4;
        real_t xd3 = C2 * wd3 + adn * yd1;
        real_t xd4 = C2 * wd4 + adn * yd2;

        real_t Fdn = 0.5 * alpha;
        if (is == 0) Fdn = -0.5 * mq;

        vt1 += Fdn * (xd1 - xd3);
        vt2 += Fdn * (xd2 - xd4);
        vt3 += Fdn * (xd3 - xd1);
        vt4 += Fdn * (xd4 - xd2);

        vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = vt4;
      }
    }
  }
}

void mult_domainwall_5din_dd_5dirdag_dirac(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c,
       int *Nsize, int *block_size, int ieo, real_t alpha){

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];

  int Nst  = Nx * Ny * Nz * Nt;

  real_t* vp_dev = (real_t*)dev_ptr(vp);
  real_t* yp_dev = (real_t*)dev_ptr(yp);
  real_t* wp_dev = (real_t*)dev_ptr(wp);

  // b, c read from __constant__ memory in the kernel (see dd_5dir_dirac).
  (void)b; (void)c;

  // Define block and grid sizes for CUDA kernel: one thread per (site, color)
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst * NVC + blockSize - 1)/ blockSize;

  mult_domainwall_5din_dd_5dirdag_dirac_dev<<<gridSize, blockSize>>>(
                   vp_dev, yp_dev, wp_dev,
                   mq, M0, Ns,
                   Nx, Ny, Nz, Nt,
                   Bx, By, Bz, Bt,
                   ieo, alpha);

  //Synchronize to ensure completion of kernel execution
  dw5din_kernel_sync();
  }

//====================================================================
__global__ void mult_domainwall_5din_dd_hopb_dirac_dev(
                            real_t *vp, real_t *up,
                            real_t *wp,
                            int Ns, 
                            int Bx, int By, int Bz, int Bt,  
                            int Nx, int Ny, int Nz, int Nt,  
                            int bc_x, int bc_y, int bc_z, int bc_t,
                            int ieo, int flag)
{
  int Nst  = Nx * Ny * Nz * Nt;
  int Nin5 = NVCD * Ns;

  // size of block
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nx/Bx;
  int NBy  = Ny/By;
  int NBz  = Nz/Bz;
  //int NBt  = Nt/Bt;
  //int Nblock = NBx * NBy * NBz * NBt;

  // Parallelize over (site, s-slice): one thread per (site2, is), like the 5D
  // hopb reference kernel. Leaving the s5 loop serial makes this kernel very
  // slow, so s5 must be parallelized here (unlike the 5d-diag kernels above).
  int gid = blockIdx.x * blockDim.x + threadIdx.x;

  if(gid < Nst * Ns){

    int is    = gid % Ns;
    int site2 = gid / Ns;

    int bsite = site2 % Bsize;
    int block = site2 /Bsize;

    int ibx = block % NBx;
    int iby = (block/NBx) % NBy;
    int ibz = (block/(NBx * NBy)) % NBz;
    int ibt = block/(NBx * NBy * NBz);

    int ieo_skip = 1 - ieo;
    int jeo = (ibx + iby + ibz + ibt) % 2;
    if(jeo == ieo_skip) return;

    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;

    int kx   = bsite % Bx;
    int ix   = kx + Bx * ibx;
    int kyzt = bsite/Bx;
    int ky   = kyzt % By;
    int iy   = ky + By * iby;
    int kzt  = kyzt/By;
    int kz   = kzt % Bz;
    int iz   = kz + Bz * ibz;
    int kt   = kzt/Bz;
    int it   = kt + Bt * ibt;
    int iyzt = iy + Ny * (iz + Nz * it);
    int izt  = iz + Nz * it;
    int ixy  = ix + Nx * iy;
    int ixyz = ix + Nx * (iy + Ny * iz);
    int site = ix + Nx * (iy + Ny * (iz + Nz * it));

    int idir;

    // thread-local gauge link buffer (replaces the racy shared sharedMemory[0],
    // which every thread overwrote). load_u rewrites it for each direction.
    real_t ut[NDF];

    {
      real_t vL[NVCD];

      if(flag == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vL[ivcd] = 0.0;
        }
      }else{
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vL[ivcd] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }
      }

      idir = 0;

      {
        int ix2 = (ix + 1) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1) bc2 = bc_x;

        load_u(ut, up, site + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_xpb(vL, ut, wt);
      }

      {
        int ix2 = (ix - 1 + Nx) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0) bc2 = bc_x;

        load_u(ut, up, nei + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
        }
        mult_wilson_xmb(vL, ut, wt);
      }

      idir = 1;

      {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc_y;

        load_u(ut, up, site + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_ypb(vL, ut, wt);
      }

      {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc_y;

        load_u(ut, up, nei + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_ymb(vL, ut, wt);
      }

      idir = 2;

      {
        int iz2 = (iz + 1) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == Nz-1) bc2 = bc_z;

        load_u(ut, up, site + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_zpb(vL, ut, wt);
      }

      {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc_z;

        load_u(ut, up, nei + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_zmb(vL, ut, wt);
      }

      idir = 3;

      {
        int it2 = (it + 1) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == Nt-1) bc2 = bc_t;

        load_u(ut, up, site + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    }
        mult_wilson_tpb_dirac(vL, ut, wt);
      }

      {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc_t;

        load_u(ut, up, nei + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	      }
        mult_wilson_tmb_dirac(vL, ut, wt);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
      }

    } // per-thread (site,is) scope
  } // if(gid < Nst*Ns)
}

void mult_domainwall_5din_dd_hopb_dirac(
                            real_t *vp, real_t *up,
                            real_t *wp,
                            int Ns, int *bc, int *Nsize,
                            int *block_size, int ieo, int flag)
{
    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Bx = block_size[0];
    int By = block_size[1];
    int Bz = block_size[2];
    int Bt = block_size[3];
    int bc_x = bc[0];
    int bc_y = bc[1];
    int bc_z = bc[2];
    int bc_t = bc[3];

    int Nst  = Nx * Ny * Nz * Nt; 
    //int Nin5 = NVCD * Ns; 

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    // one thread per (site, s-slice); gauge link is now a thread-local array,
    // so no dynamic shared memory is needed.
    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst * Ns + blockSize - 1)/ blockSize;

    // Launch the CUDA kernel
    mult_domainwall_5din_dd_hopb_dirac_dev<<<gridSize, blockSize>>>(
           vp_dev, up_dev, wp_dev,
           Ns,
           Bx, By, Bz, Bt,
           Nx, Ny, Nz, Nt,
           bc_x, bc_y, bc_z, bc_t,
           ieo, flag);

    //Synchronize to ensure completion of kernel execution
    dw5din_kernel_sync();
};

#endif
//============================================================END=====
