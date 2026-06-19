/*!
       @file    mult_Doainwall_5din_cuda2-inc.h
       @brief
       @author  Wei-Lun Chen (wlchen)
                $LastChangedBy: matufuru $
       @date    $LastChangedDate: 2024-03-14 20:51:53 #$
       @version $LastChangedRevision: 2581 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC2_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC2_INCLUDED

//====================================================================
__global__ void mult_domainwall_5din_xpb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int idir = 0;
    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((ix < Nx-1) || (do_comm == 0)) {
                int ix2 = (ix + 1) % Nx;
                int nei = ix2 + Nx * iyzt;
                real_t bc2 = (ix == Nx - 1) ? bc : 1.0;
                real_t* ut = &sharedMemory[0];
                load_u(ut, up, site + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_xpb(vL, ut, wt);
            }
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_xpb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[0];
    int bc_dev      = bc[0];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_xpb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ void mult_domainwall_5din_xmb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int idir = 0;

    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((ix > 0) || (do_comm == 0)) {

                int ix2 = (ix - 1 + Nx) % Nx;
                int nei = ix2 + Nx * iyzt;

                real_t bc2 = (ix == 0) ? bc : 1.0;
                real_t* ut = &sharedMemory[0];
                load_u(ut, up, nei + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_xmb(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_xmb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[0];
    int bc_dev      = bc[0];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_xmb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ void mult_domainwall_5din_ypb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int idir = 1;

    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int iy   = iyzt % Ny;
        int izt  = iyzt / Ny;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((iy < Ny-1) || (do_comm == 0)) {
                int iy2 = (iy + 1) % Ny;
                int nei = ix + Nx * (iy2 + Ny * izt);
                real_t bc2 = (iy == Ny-1) ? bc : 1.0;

                real_t* ut = &sharedMemory[0];
                load_u(ut, up, site + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_ypb(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_ypb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[1];
    int bc_dev      = bc[1];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_ypb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}

//====================================================================
__global__ void mult_domainwall_5din_ymb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx  * Ny;
    int idir = 1;

    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ix   = site % Nx;
        int iyzt = site / Nx;
        int iy   = iyzt % Ny;
        int izt  = iyzt / Ny;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((iy > 0) || (do_comm == 0)) {
                int iy2 = (iy - 1 + Ny) % Ny;
                int nei = ix + Nx * (iy2 + Ny * izt);
                real_t bc2 = (iy == 0) ? bc : 1.0;

                real_t* ut = &sharedMemory[0];
                load_u(ut, up, nei + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_ymb(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_ymb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[1];
    int bc_dev      = bc[1];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_ymb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}
//====================================================================
__global__ void mult_domainwall_5din_zpb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx  * Ny;
    int idir = 2;
    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ixy  = site % Nxy;
        int izt  = site / Nxy;
        int iz   = izt % Nz;
        int it   = izt / Nz;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((iz < Nz-1) || (do_comm == 0)) {
                int iz2 = (iz + 1) % Nz;
                int nei = ixy + Nxy * (iz2 + Nz * it);
                real_t bc2 = (iz == Nz-1) ? bc : 1.0;

                real_t* ut = &sharedMemory[0];
                load_u(ut, up, site + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_zpb(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_zpb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[2];
    int bc_dev      = bc[2];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_zpb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}
//====================================================================
__global__ void mult_domainwall_5din_zmb_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx  * Ny;
    int idir = 2;
    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int ixy  = site % Nxy;
        int izt  = site / Nxy;
        int iz   = izt % Nz;
        int it   = izt / Nz;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((iz > 0) || (do_comm == 0)) {
                int iz2 = (iz - 1 + Nz) % Nz;
                int nei = ixy + Nxy * (iz2 + Nz * it);

                real_t bc2 = (iz == 0) ? bc : 1.0;
                real_t* ut = &sharedMemory[0];
                load_u(ut, up, nei + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_zmb(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_zmb(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[2];
    int bc_dev      = bc[2];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_zmb_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}
//====================================================================
__global__ void mult_domainwall_5din_tpb_dirac_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;
    int idir = 3;

    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int izt  = site  / Nxy;
        int it   = izt   / Nz;
        int ixyz = site  % Nxyz;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((it < Nt-1) || (do_comm == 0)) {
                int it2 = (it + 1) % Nt;
                int nei = ixyz + Nxyz * it2;

                real_t bc2 = (it == Nt-1) ? bc : 1.0;
                real_t* ut = &sharedMemory[0];
                load_u(ut, up, site + Nst * idir);

                real_t wt[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	            }
                mult_wilson_tpb_dirac(vL, ut, wt);
            }

            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vL[ivcd];
            }
        } // is loop
    }
}

void mult_domainwall_5din_tpb_dirac(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[3];
    int bc_dev      = bc[3];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_tpb_dirac_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}
//====================================================================
__global__ void mult_domainwall_5din_tmb_dirac_dev(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int bc, 
     int Nx, int Ny, int Nz, int Nt, 
     int do_comm, int flag)
{
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    int Nst  = Nx * Ny * Nz * Nt;
    int Nin5 = NVCD * Ns;
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;
    int idir = 3;

    extern __shared__ real_t sharedMemory[];

    if (site < Nst) {

        int izt  = site / Nxy;
        int it   = izt  / Nz;
        int ixyz = site % Nxyz;

        for(int is = 0; is < Ns; ++is){

            real_t vL[NVCD];

            for(int ivcd = 0; ivcd < NVCD; ++ivcd) {
                vL[ivcd] = (flag == 0) ? 0.0 : vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }

            if ((it > 0) || (do_comm == 0)) {
                int it2 = (it - 1 + Nt) % Nt;
                int nei = ixyz + Nxyz * it2;

                real_t bc2 = (it == 0) ? bc : 1.0;
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
    }
}

void mult_domainwall_5din_tmb_dirac(
     real_t *vp, real_t *up, real_t *wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{

    int Nx = Nsize[0];
    int Ny = Nsize[1];
    int Nz = Nsize[2];
    int Nt = Nsize[3];
    int Nst = Nx * Ny * Nz * Nt; 

    int do_comm_dev = do_comm[3];
    int bc_dev      = bc[3];

    real_t* vp_dev = (real_t*)dev_ptr(vp);
    real_t* up_dev = (real_t*)dev_ptr(up);
    real_t* wp_dev = (real_t*)dev_ptr(wp);

    int blockSize  = VECTOR_LENGTH;
    int gridSize   = (Nst + blockSize - 1)/ blockSize;
    int sharedsize = NDF * blockSize * sizeof(real_t);

    // Launch the CUDA kernel
    mult_domainwall_5din_tpb_dirac_dev<<<gridSize, blockSize, sharedsize>>>(
    vp_dev, up_dev, wp_dev, Ns, bc_dev, 
    Nx, Ny, Nz, Nt, do_comm_dev, flag);

    // Synchronize to ensure completion of kernel execution
    CHECK(cudaDeviceSynchronize());
}

__global__  void mult_domainwall_5din_xp1_dev(
        real_t * buf_xp,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int iyzt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nyzt         = Ny   * Nz  * Nt     ;
        int Nin5         = NVCD * Ns           ;
        int Nin5bd       = NVC  * ND2 * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int idir         = 0                   ;

        if ( iyzt < Nyzt ){

            int ix   = 0;
            int site = ix + Nx * iyzt;
            real_t bc2 = bc;

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_xp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)] = bc2 * vt[ivcd];
                }
            } // is loop
        }
    }

__global__  void mult_domainwall_5din_xm1_dev(
        real_t *buf_xm,
        real_t *up, real_t *wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int iyzt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nyzt         = Ny   * Nz  * Nt     ;
        int Nin5         = NVCD * Ns           ;
        int Nin5bd       = NVC  * ND2 * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int idir         = 0                   ;
        extern __shared__ real_t sharedMemory[];

        if ( iyzt < Nyzt ){

            int ix   = Nx - 1;
            int site = ix + Nx * iyzt;

            real_t* ut = &sharedMemory[0];
            load_u(ut, up, site + Nst * idir);

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_xm1(vt, ut, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)] = vt[ivcd];
                }
            }//is loop
        }
    }

__global__  void mult_domainwall_5din_yp1_dev(
        real_t * buf_yp,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixzt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxzt         = Nx * Nz * Nt;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir = 1;

        if ( ixzt < Nxzt ){

            int iy   = 0;
            int ix   = ixzt % Nx;
            int izt  = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);
            real_t bc2 = bc;

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_yp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = bc2 * vt[ivcd];
                }
            }// is loop
        }
    }

__global__  void mult_domainwall_5din_ym1_dev(
        real_t * buf_ym,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixzt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxzt         = Nx * Nz * Nt;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir = 1;
        extern __shared__ real_t sharedMemory[];

        if ( ixzt < Nxzt ){
            
            int iy   = Ny - 1;    
            int ix   = ixzt % Nx;
            int izt  = ixzt / Nx;
            int site = ix + Nx * (iy + Ny * izt);

            real_t bc2 = bc;
            real_t* ut = &sharedMemory[0];
            load_u(ut, up , site + Nst * idir);

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_ym1(vt, ut, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = vt[ivcd];
                }
            } // is loop
        }
    }


__global__  void mult_domainwall_5din_zp1_dev(
        real_t * buf_zp,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixyt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nxyt         = Nx * Ny * Nt;
        int Nxy          = Nx * Ny;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir = 2;

        if ( ixyt < Nxyt ){

            int iz    = 0;
            int ixy   = ixyt % Nxy;
            int it    = ixyt / Nxy;
            int site  = ixy  + Nxy * (iz + Nz * it);
            real_t bc2 = bc;

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_zp1(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = bc2 * vt[ivcd];
                }
            }// is loop
        }
    }

__global__  void mult_domainwall_5din_zm1_dev(
        real_t * buf_zm,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixyt         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx * Ny;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nxyt         = Nx * Ny * Nt;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir = 2;
        extern __shared__ real_t sharedMemory[];

        if ( ixyt < Nxyt ){

            int ixy  = ixyt % Nxy;
            int it   = ixyt / Nxy;
            int iz   = Nz - 1;
            int site = ixy + Nxy * (iz + Nz * it);
            real_t* ut = &sharedMemory[0];
            load_u(ut, up , site + Nst * idir);

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_zm1(vt, ut, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = vt[ivcd];
                }
            } // is loop
        }
    }

__global__  void mult_domainwall_5din_tp1_dirac_dev(
        real_t * buf_tp,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixyz         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxyz         = Nx * Ny * Nz;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir   = 3;

        if ( ixyz < Nxyz ){

            int it     = 0;
            int site   = ixyz + Nxyz * it;
            real_t bc2 = bc;

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_tp1_dirac(vt, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = bc2 * vt[ivcd];
                }
            } // is loop
        }
    }

__global__  void mult_domainwall_5din_tm1_dirac_dev(
        real_t * buf_tm,
        real_t * up, real_t * wp,
        int Ns, int bc, int Nx, int Ny, int Nz, int Nt){

        int ixyz         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nst          = Nx * Ny * Nz * Nt;
        int Nxyz         = Nx * Ny * Nz;
        int Nin5         = NVCD * Ns;
        int Nin5bd       = NVC  * ND2 * Ns;
        int idir         = 3; 
        extern __shared__ real_t sharedMemory[];

        if ( ixyz < Nxyz ){

            int it    = Nt - 1;
            int site  = ixyz + Nxyz * it;
                    
            real_t* ut = &sharedMemory[0];
            load_u(ut, up , site + Nst *idir);

            for (int is = 0; is < Ns; ++is) {
                real_t wt[NVCD], vt[NVC * ND2];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
                }
                mult_wilson_tm1_dirac(vt, ut, wt);
                for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
                    buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = vt[ivcd];
                }
            }// is loop
        }
    }

void mult_domainwall_5din_xp1(
        real_t *buf1_xp, real_t *up, real_t *wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstx   = Ny * Nz * Nt;
        int bc_x   = bc[0];

        real_t * buf1_xp_dev = (real_t*)dev_ptr(buf1_xp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);

        int blockSize  =  VECTOR_LENGTH;
        int gridSize   =  (Nstx + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_xp1_dev<<<gridSize, blockSize>>>(
        buf1_xp_dev,
        up_dev, wp_dev,
        Ns, bc_x, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_xm1(
        real_t * buf1_xm, real_t * up, real_t * wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstx   = Ny * Nz * Nt;
        int bc_x   = bc[0];

        real_t * buf1_xm_dev = (real_t*)dev_ptr(buf1_xm);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);

        int blockSize  =  VECTOR_LENGTH;
        int gridSize   =  (Nstx + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_xm1_dev<<<gridSize, blockSize, sharedsize>>>(
        buf1_xm_dev,
        up_dev, wp_dev,
        Ns, bc_x, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_yp1(
        real_t *buf1_yp, real_t *up, real_t *wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nsty   = Nz * Nz * Nt;
        int bc_y   = bc[1];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_yp_dev = (real_t*)dev_ptr(buf1_yp);

        int blockSize  =  VECTOR_LENGTH;
        int gridSize   =  (Nsty + blockSize - 1)/ blockSize;

        mult_domainwall_5din_yp1_dev<<<gridSize, blockSize>>>(
        buf1_yp_dev,
        up_dev, wp_dev,
        Ns, bc_y, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_ym1(
        real_t * buf1_ym, real_t * up, real_t * wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nsty   = Nx * Nz * Nt;
        int bc_y   = bc[1];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_ym_dev = (real_t*)dev_ptr(buf1_ym);

        int blockSize =  VECTOR_LENGTH;
        int gridSize  =  (Nsty + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_ym1_dev<<<gridSize, blockSize, sharedsize>>>(
        buf1_ym_dev,
        up_dev, wp_dev,
        Ns, bc_y, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_zp1(
        real_t *buf1_zp, real_t *up, real_t *wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstz   = Nx * Ny * Nt;
        int bc_z   = bc[2];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_zp_dev = (real_t*)dev_ptr(buf1_zp);

        int blockSize =  VECTOR_LENGTH;
        int gridSize  =  (Nstz + blockSize - 1)/ blockSize;

        mult_domainwall_5din_zp1_dev<<<gridSize, blockSize>>>(
        buf1_zp_dev,
        up_dev, wp_dev,
        Ns, bc_z, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_zm1(
        real_t * buf1_zm, real_t * up, real_t * wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstz   = Nx * Ny * Nt;
        int bc_z   = bc[2];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_zm_dev = (real_t*)dev_ptr(buf1_zm);

        int blockSize =  VECTOR_LENGTH;
        int gridSize  =  (Nstz + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_zm1_dev<<<gridSize, blockSize, sharedsize>>>(
        buf1_zm_dev,
        up_dev, wp_dev,
        Ns, bc_z, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_tp1_dirac(
        real_t *buf1_tp, real_t *up, real_t *wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstt   = Nx * Ny * Nz;
        int bc_t   = bc[3];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_tp_dev = (real_t*)dev_ptr(buf1_tp);

        int blockSize =  VECTOR_LENGTH;
        int gridSize  =  (Nstt + blockSize - 1)/ blockSize;

        mult_domainwall_5din_tp1_dirac_dev<<<gridSize, blockSize>>>(
        buf1_tp_dev,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

void mult_domainwall_5din_tm1_dirac(
        real_t * buf1_tm, real_t * up, real_t * wp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nstt   = Nx * Ny * Nz;
        int bc_t   = bc[3];

        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * wp_dev = (real_t*)dev_ptr(wp);
        real_t * buf1_tm_dev = (real_t*)dev_ptr(buf1_tm);

        int blockSize  =  VECTOR_LENGTH;
        int gridSize   =  (Nstt + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_tm1_dirac_dev<<<gridSize, blockSize, sharedsize>>>(
        buf1_tm_dev,
        up_dev, wp_dev,
        Ns, bc_t, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
};

__global__ void mult_domainwall_5din_xp2_dev(
    real_t *vp, real_t *up,
    real_t *buf_xp, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;
        int idir         = 0;
        extern __shared__ real_t sharedMemory[];

        if (site < Nst) {

            int ix   = site % Nx;
            int iyzt = site / Nx;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                idir = 0;

                if (ix == Nx-1) {
                    real_t* ut = &sharedMemory[0];
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xp2(vL, ut, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }


__global__ void mult_domainwall_5din_xm2_dev(
    real_t *vp, real_t *up,
    real_t *buf_xm, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;
        int idir         = 0;

        if (site < Nst) {

            int ix   = site % Nx;
            int iyzt = site / Nx;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                idir = 0;

                if (ix == 0) {
                    real_t bc2 = bc;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
                    }
                    mult_wilson_xm2(vL, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }


__global__ void mult_domainwall_5din_yp2_dev(
    real_t *vp, real_t *up,
    real_t *buf_yp, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;
        extern __shared__ real_t sharedMemory[];

        if (site < Nst) {

            int ix   = site % Nx;
            int iyzt = site / Nx;
            int ixy  = site % Nxy;
            int iy   = iyzt % Ny;
            int izt  = site / Nxy;
            int idir = 1;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                int ixzt = ix + Nx * izt;

                if (iy == Ny-1) {
                    real_t* ut = &sharedMemory[0];
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC * ND2 * is), ixzt)];
                    }
                    mult_wilson_yp2(vL, ut, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            } // is loop
        }
    }

__global__ void mult_domainwall_5din_ym2_dev(
    real_t *vp, real_t *up,
    real_t *buf_ym, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;

        if (site < Nst) {

            int ix   = site % Nx;
            int iyzt = site / Nx;
            int ixy  = site % Nxy;
            int iy   = iyzt % Ny;
            int izt  = site / Nxy;
            int idir = 1;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }
                int ixzt = ix + Nx * izt;

                if (iy == 0) {
                    real_t bc2 = bc;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
                    }
                    mult_wilson_ym2(vL, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }

__global__ void mult_domainwall_5din_zp2_dev(
    real_t *vp, real_t *up,
    real_t *buf_zp, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;
        extern __shared__ real_t sharedMemory[];

        if (site < Nst) {

            int ixy  = site % Nxy;
            int izt  = site / Nxy;
            int iz   = izt  % Nz;
            int it   = izt  / Nz;
            int idir = 2;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];
                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                int ixyt = ixy + Nxy * it;

                if (iz == Nz-1) {
                    real_t* ut = &sharedMemory[0];
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zp2(vL, ut, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }

__global__ void mult_domainwall_5din_zm2_dev(
    real_t *vp, real_t *up,
    real_t *buf_zm, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;

        if (site < Nst) {

            int ixy  = site % Nxy;
            int izt  = site / Nxy;
            int iz   = izt  % Nz;
            int it   = izt  / Nz;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                int ixyt = ixy + Nxy * it;

                if (iz == 0) {
                    real_t bc2 = bc;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
                    }
                    mult_wilson_zm2(vL, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }

__global__ void mult_domainwall_5din_tp2_dirac_dev(
    real_t *vp, real_t *up,
    real_t *buf_tp, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;
        extern __shared__ real_t sharedMemory[];

        if (site < Nst) {

            int izt  = site / Nxy;
            int it   = izt  / Nz;
            int ixyz = site % Nxyz;
            int idir = 3;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }
    
                if (it == Nt-1) {
                    real_t* ut = &sharedMemory[0];
                    load_u(ut, up, site + Nst * idir);
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tp2_dirac(vL, ut, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }


__global__ void mult_domainwall_5din_tm2_dirac_dev(
    real_t *vp, real_t *up,
    real_t *buf_tm, 
    int Ns, int bc, 
    int Nx, int Ny, int Nz, int Nt) {

        int site         = blockIdx.x * blockDim.x + threadIdx.x;
        int Nxy          = Nx   * Ny           ;
        int Nxyz         = Nx   * Ny  * Nz     ;
        int Nin5bd       = (NVCD / 2) * Ns     ;
        int Nst          = Nx   * Ny  * Nz * Nt;
        int Nin5         = NVCD * Ns           ;

        if (site < Nst) {

            int izt  = site / Nxy;
            int it   = izt  / Nz;
            int ixyz = site % Nxyz;
            int idir = 3;

            for(int is = 0; is < Ns; ++is){

                real_t vL[NVCD];

                for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                    vL[ivcd] = 0.0;
                }

                if (it == 0) {
                    real_t bc2 = bc;
                    real_t wt[NVC * ND2];
                    for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
                        wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
                    }
                    mult_wilson_tm2_dirac(vL, wt);
                    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
                        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	                }
                }
            }
        }
    }


void mult_domainwall_5din_xp2(
        real_t *vp     , real_t *up,
        real_t *buf2_xp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 0;
        int bc_x   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_xp_dev = (real_t*)dev_ptr(buf2_xp);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_xp2_dev<<<gridSize, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,
            buf2_xp_dev,
            Ns,  bc_x, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}
void mult_domainwall_5din_xm2(
        real_t *vp     , real_t *up,
        real_t *buf2_xm,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 0;
        int bc_x   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_xm_dev = (real_t*)dev_ptr(buf2_xm);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;

        mult_domainwall_5din_xp2_dev<<<gridSize, VECTOR_LENGTH>>>(
            vp_dev,  up_dev,
            buf2_xm_dev,
            Ns,  bc_x, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}
void mult_domainwall_5din_yp2(
        real_t *vp     , real_t *up,
        real_t *buf2_yp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 1;
        int bc_y = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_yp_dev = (real_t*)dev_ptr(buf2_yp);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_yp2_dev<<<gridSize, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,
            buf2_yp_dev,
            Ns,  bc_y, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_ym2(
        real_t *vp     , real_t *up,
        real_t *buf2_ym,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 1;
        int bc_y = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_ym_dev = (real_t*)dev_ptr(buf2_ym);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;

        mult_domainwall_5din_ym2_dev<<<gridSize, VECTOR_LENGTH>>>(
            vp_dev,  up_dev,
            buf2_ym_dev,
            Ns,  bc_y, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_zp2(
        real_t *vp     , real_t *up,
        real_t *buf2_zp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 2;
        int bc_z   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_zp_dev = (real_t*)dev_ptr(buf2_zp);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_zp2_dev<<<gridSize, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev, 
            buf2_zp_dev,
            Ns,  bc_z, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_zm2(
        real_t *vp, real_t *up, real_t *buf2_zm,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 2; 
        int bc_z   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_zm_dev = (real_t*)dev_ptr(buf2_zm);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;

        mult_domainwall_5din_zm2_dev<<<gridSize, VECTOR_LENGTH>>>(
            vp_dev,  up_dev, 
            buf2_zm_dev,
            Ns,  bc_z, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_tp2_dirac(
        real_t *vp     , real_t *up,
        real_t *buf2_tp,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 3;
        int bc_t   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_tp_dev = (real_t*)dev_ptr(buf2_tp);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;
        int sharedsize = NDF * blockSize * sizeof(real_t);

        mult_domainwall_5din_tp2_dirac_dev<<<gridSize, VECTOR_LENGTH, sharedsize>>>(
            vp_dev,  up_dev,
            buf2_tp_dev,
            Ns,  bc_t, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

void mult_domainwall_5din_tm2_dirac(
        real_t *vp , real_t *up,
        real_t *buf2_tm,
        int Ns, int *bc, int *Nsize){

        int Nx     = Nsize[0];
        int Ny     = Nsize[1];
        int Nz     = Nsize[2];
        int Nt     = Nsize[3];
        int Nst    = Nx * Ny * Nz * Nt;
        int idir   = 3;
        int bc_t   = bc[idir];

        real_t * vp_dev = (real_t*)dev_ptr(vp);
        real_t * up_dev = (real_t*)dev_ptr(up);
        real_t * buf2_tm_dev = (real_t*)dev_ptr(buf2_tm);

        int blockSize  = VECTOR_LENGTH;
        int gridSize   = (Nst + blockSize - 1)/ blockSize;

        mult_domainwall_5din_tm2_dirac_dev<<<gridSize, VECTOR_LENGTH>>>(
            vp_dev,  up_dev,
            buf2_tm_dev,
            Ns,  bc_t, Nx, Ny, Nz, Nt);

        CHECK(cudaDeviceSynchronize());
}

#endif
//============================================================END=====
