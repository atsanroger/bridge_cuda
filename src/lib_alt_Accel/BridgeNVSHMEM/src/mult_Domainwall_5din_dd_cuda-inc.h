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
                   real_t *b_dev, real_t *c_dev,
                   int Nx, int Ny, int Nz, int Nt,
                   int Bx, int By, int Bz, int Bt,
                   int ieo)
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

  int site2  = blockIdx.x * blockDim.x + threadIdx.x;

  extern __shared__ real_t shared_mem[];
  real_t* b = shared_mem;
  real_t* c = &shared_mem[Ns]; 

  int tid = threadIdx.x;
  if(tid < Ns) {
    b[tid] = b_dev[tid];
    c[tid] = c_dev[tid];
  }

  if (site2 < Nst) {

    int bsite = site2 % Bsize;
    int block = site2 / Bsize;

    int ibx = block % NBx;
    int iby = (block/NBx) % NBy;
    int ibz = (block/(NBx * NBy)) % NBz;
    int ibt = block/(NBx * NBy * NBz);

    int ieo_skip = 1 - ieo;
    int jeo = (ibx + iby + ibz + ibt) % 2;
    if(jeo != ieo_skip){

      int Nxy  = Nx  * Ny;
      int Nxyz = Nxy * Nz;

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
      int iyzt = iy + Ny * (iz + Nz * it);
      int izt  = iz + Nz * it;
      int ixy  = ix + Nx * iy;
      int ixyz = ix + Nx * (iy + Ny * iz);
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      for (int is = 0; is < Ns; ++is) {

        real_t vt[NVCD], wt[NVCD];

        int is_up = (is+1) % Ns;
        real_t Fup   = 0.5;
        if (is == Ns-1) Fup = -0.5 * mq;

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd]  = wp[IDX2(Nin5, (ivcd + NVCD*is_up), site)];
        }
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vt[ID1 + ivc]  = Fup * (wt[ID1 + ivc] - wt[ID3 + ivc]);
          vt[ID2 + ivc]  = Fup * (wt[ID2 + ivc] - wt[ID4 + ivc]);
          vt[ID3 + ivc]  = Fup * (wt[ID3 + ivc] - wt[ID1 + ivc]);
          vt[ID4 + ivc]  = Fup * (wt[ID4 + ivc] - wt[ID2 + ivc]);
        }

        int is_dn = (is-1 + Ns) % Ns;
        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is_dn), site)];
        }
        for (int ivc = 0; ivc < NVC; ++ivc) {
          vt[ID1 + ivc] += Fdn * (wt[ID1 + ivc] + wt[ID3 + ivc]);
          vt[ID2 + ivc] += Fdn * (wt[ID2 + ivc] + wt[ID4 + ivc]);
          vt[ID3 + ivc] += Fdn * (wt[ID3 + ivc] + wt[ID1 + ivc]);
          vt[ID4 + ivc] += Fdn * (wt[ID4 + ivc] + wt[ID2 + ivc]);
        }

        real_t B_is = b[is] * (4.0 - M0) + 1.0;
        real_t C_is = c[is] * (4.0 - M0) - 1.0;

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
        }

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          vp[IDX2(Nin5, (ivcd + NVCD * is), site)]
                = B_is * wt[ivcd] + C_is * vt[ivcd];
          yp[IDX2(Nin5, (ivcd + NVCD * is), site)]
                = -0.5 * (b[is] * wt[ivcd] + c[is] * vt[ivcd]);
        }
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
                   int *Nsize, int *block_size, int ieo)
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

  real_t* b_dev;
  real_t* c_dev; 

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&b_dev , nsize));
  CHECK(cudaMalloc((void **)&c_dev , nsize));

  CHECK(cudaMemcpy(b_dev, b, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(c_dev, c, nsize, cudaMemcpyHostToDevice));

  // Define block and grid sizes for CUDA kernel
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst + blockSize - 1)/ blockSize;
  int shared_mem_size = 2 * nsize;

  mult_domainwall_5din_dd_5dir_dirac_dev<<<gridSize, blockSize,
                                           shared_mem_size>>>(
                   vp_dev, yp_dev, wp_dev,
                   mq, M0, Ns, 
                   b_dev, c_dev,
                   Nx, Ny, Nz, Nt,
                   Bx, By, Bz, Bt,
                   ieo);
              
  //Synchronize to ensure completion of kernel execution
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(b_dev));
  CHECK(cudaFree(c_dev));
}

//====================================================================
__global__ void mult_domainwall_5din_dd_5dirdag_dirac_dev(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns, 
       real_t *b_dev, real_t *c_dev,
       int Nx, int Ny, int Nz, int Nt,
       int Bx, int By, int Bz, int Bt,
       int ieo)
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

  extern __shared__ real_t shared_mem[];
  real_t* b = shared_mem;
  real_t* c = &shared_mem[Ns]; 

  int tid = threadIdx.x;
  if(tid < Ns) {
    b[tid] = b_dev[tid];
    c[tid] = c_dev[tid];
  }

  int site2  = blockIdx.x * blockDim.x + threadIdx.x;

  if( site2 < Nst ){

    int bsite = site2 % Bsize;
    int block = site2/Bsize;
    int ibx = block % NBx;
    int iby = (block/NBx) % NBy;
    int ibz = (block/(NBx * NBy)) % NBz;
    int ibt = block/(NBx * NBy * NBz);

    int ieo_skip = 1 - ieo;
    int jeo = (ibx + iby + ibz + ibt) % 2;
    if(jeo != ieo_skip){

      int Nxy  = Nx  * Ny;
      int Nxyz = Nxy * Nz;

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
      int iyzt = iy + Ny * (iz + Nz * it);
      int izt  = iz + Nz * it;
      int ixy  = ix + Nx * iy;
      int ixyz = ix + Nx * (iy + Ny * iz);
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      for (int is = 0; is < Ns; ++is) {
        real_t vt[NVCD], xt[NVCD], wt[NVCD], yt[NVCD];
        real_t B1 = b[is] * (4.0 - M0) + 1.0;
        real_t a1 = -0.5 * b[is];

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is), site)];
        }

        for (int ivc = 0; ivc < NVC; ++ivc) {
          vt[ID1 + ivc] = B1 * wt[ID1 + ivc] + a1 * yt[ID3 + ivc];
          vt[ID2 + ivc] = B1 * wt[ID2 + ivc] + a1 * yt[ID4 + ivc];
          vt[ID3 + ivc] = B1 * wt[ID3 + ivc] + a1 * yt[ID1 + ivc];
          vt[ID4 + ivc] = B1 * wt[ID4 + ivc] + a1 * yt[ID2 + ivc];
        }

        int is_up = (is+1) % Ns;
        real_t C1 = c[is_up] * (4.0 - M0) - 1.0;
        real_t aup = -0.5 * c[is_up];

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
          yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
        }

        for (int ivc = 0; ivc < NVC; ++ivc) {
          xt[ID1 + ivc] = C1 * wt[ID1 + ivc] + aup * yt[ID3 + ivc];
          xt[ID2 + ivc] = C1 * wt[ID2 + ivc] + aup * yt[ID4 + ivc];
          xt[ID3 + ivc] = C1 * wt[ID3 + ivc] + aup * yt[ID1 + ivc];
          xt[ID4 + ivc] = C1 * wt[ID4 + ivc] + aup * yt[ID2 + ivc];
        }

        real_t Fup = 0.5;
        if (is == Ns-1) Fup = -0.5 * mq;

        for (int ivc = 0; ivc < NVC; ++ivc) {
          vt[ID1 + ivc] += Fup * (xt[ID1 + ivc] + xt[ID3 + ivc]);
          vt[ID2 + ivc] += Fup * (xt[ID2 + ivc] + xt[ID4 + ivc]);
          vt[ID3 + ivc] += Fup * (xt[ID3 + ivc] + xt[ID1 + ivc]);
          vt[ID4 + ivc] += Fup * (xt[ID4 + ivc] + xt[ID2 + ivc]);
        }

        int is_dn = (is-1 + Ns) % Ns;
        real_t C2 = c[is_dn] * (4.0 - M0) - 1.0;
        real_t adn = -0.5 * c[is_dn];

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
          yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
        }

        for (int ivc = 0; ivc < NVC; ++ivc) {
          xt[ID1 + ivc] = C2 * wt[ID1 + ivc] + adn * yt[ID3 + ivc];
          xt[ID2 + ivc] = C2 * wt[ID2 + ivc] + adn * yt[ID4 + ivc];
          xt[ID3 + ivc] = C2 * wt[ID3 + ivc] + adn * yt[ID1 + ivc];
          xt[ID4 + ivc] = C2 * wt[ID4 + ivc] + adn * yt[ID2 + ivc];
        }

        real_t Fdn   = 0.5;
        if (is == 0) Fdn = -0.5 * mq;

        for (int ivc = 0; ivc < NVC; ++ivc) {
          vt[ID1 + ivc] += Fdn * (xt[ID1 + ivc] - xt[ID3 + ivc]);
          vt[ID2 + ivc] += Fdn * (xt[ID2 + ivc] - xt[ID4 + ivc]);
          vt[ID3 + ivc] += Fdn * (xt[ID3 + ivc] - xt[ID1 + ivc]);
          vt[ID4 + ivc] += Fdn * (xt[ID4 + ivc] - xt[ID2 + ivc]);
        }

        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
        }
      }
    }
  }
}

void mult_domainwall_5din_dd_5dirdag_dirac(
       real_t *vp, real_t *yp, real_t *wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c,
       int *Nsize, int *block_size, int ieo){

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

  real_t* b_dev;
  real_t* c_dev; 

  size_t nsize = Ns * sizeof(real_t);

  CHECK(cudaMalloc((void **)&b_dev , nsize));
  CHECK(cudaMalloc((void **)&c_dev , nsize));

  CHECK(cudaMemcpy(b_dev, b, nsize, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(c_dev, c, nsize, cudaMemcpyHostToDevice));

  // Define block and grid sizes for CUDA kernel
  int blockSize = VECTOR_LENGTH;
  int gridSize  = (Nst + blockSize - 1)/ blockSize;
  int shared_mem_size = 2 * nsize;

  mult_domainwall_5din_dd_5dirdag_dirac_dev<<<gridSize, blockSize,
                                              shared_mem_size>>>(
                   vp_dev, yp_dev, wp_dev,
                   mq, M0, Ns, 
                   b_dev, c_dev,
                   Nx, Ny, Nz, Nt,
                   Bx, By, Bz, Bt,
                   ieo);
              
  //Synchronize to ensure completion of kernel execution
  CHECK(cudaDeviceSynchronize());
  CHECK(cudaFree(b_dev));
  CHECK(cudaFree(c_dev));
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

  int site2  = blockIdx.x * blockDim.x + threadIdx.x;

  extern __shared__ real_t sharedMemory[];

  if(site2 < Nst){

    int bsite = site2 % Bsize;
    int block = site2/Bsize;

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

    for(int is = 0; is < Ns; ++is){

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

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
        load_u(ut, up, nei + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	    
        mult_wilson_xmb(vL, ut, wt);
      }

      idir = 1;

      {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc_y;

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
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

        real_t* ut = &sharedMemory[0];
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

    } // is loop
   } // site loop
  }
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

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_dd_hopb_dirac_dev<<<gridSize, blockSize,
                                             sharedsize>>>(
           vp_dev, up_dev, wp_dev, 
           Ns, 
           Bx, By, Bz, Bt,
           Nx, Ny, Nz, Nt, 
           bc_x, bc_y, bc_z, bc_t,
           ieo, flag);

    //Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
};

#endif
//============================================================END=====
