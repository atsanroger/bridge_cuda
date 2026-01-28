/*!
      @file    mult_Doainwall_5din_openacc2-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2575 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC2_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC2_INCLUDED

//====================================================================
void mult_domainwall_5din_xpb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;

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

      if ((ix < Nx-1) || (do_comm[idir] == 0)) {
        int ix2 = (ix + 1) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1) bc2 = bc[0];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_xmb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {

    int ix   = site % Nx;
    int iyzt = site / Nx;

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

      if ((ix > 0) || (do_comm[idir] == 0)) {
        int ix2 = (ix - 1 + Nx) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0) bc2 = bc[0];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_ypb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int iy   = iyzt % Ny;
    int izt  = iyzt / Ny;

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

      int idir = 1;

      if ((iy < Ny-1) || (do_comm[idir] == 0)) {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc[idir];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_ymb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;

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

      idir = 1;

      if ((iy > 0) || (do_comm[idir] == 0)) {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc[1];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_zpb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;

    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

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

      idir = 2;

      if ((iz < Nz-1) || (do_comm[idir] == 0)) {
        int iz2 = (iz + 1) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == Nz-1) bc2 = bc[2];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_zmb(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;

    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

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

      idir = 2;

      if ((iz > 0) || (do_comm[idir] == 0)) {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc[2];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_tpb_dirac(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;

    int izt  = site / Nxy;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

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

      idir = 3;

      if ((it < Nt-1) || (do_comm[idir] == 0)) {
        int it2 = (it + 1) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == Nt-1) bc2 = bc[3];

        real_t ut[NDF];
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

  } // site loop
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_tmb_dirac(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;

    int izt  = site / Nxy;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

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

      idir = 3;

      if ((it > 0) || (do_comm[idir] == 0)) {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc[3];

        real_t ut[NDF];
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
 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_xp1(real_t *RESTRICT buf_xp,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size   = Nin5 * Nst;
  int size_u = NDF * Nst * NDIM;

  int size_bx = Nin5bd * (Nst/Nx);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xp[0:size_bx])				\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int Nyzt = Ny * Nz * Nt;

#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int iyzt = 0; iyzt < Nyzt; ++iyzt) {
    int ix = 0;
    int site = ix + Nx * iyzt;
    real_t bc2 = bc[0];
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_xp1(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)] = bc2 * vt[ivcd];
      }
    }
  }

 }
}

//====================================================================
void mult_domainwall_5din_xm1(real_t *RESTRICT buf_xm,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_bx = Nin5bd * (Nst/Nx);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xm[0:size_bx])				\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int idir = 0;
  int Nyzt = Ny * Nz * Nt;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int iyzt = 0; iyzt < Nyzt; ++iyzt) {
    int ix = Nx-1;
    int site = ix + Nx * iyzt;
    real_t ut[NDF];
    load_u(ut, up, site + Nst*idir);
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_xm1(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)] = vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_yp1(real_t *RESTRICT buf_yp,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_by = Nin5bd * (Nst/Ny);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_yp[0:size_by])				\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int Nxzt = Nx * Nz * Nt;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixzt = 0; ixzt < Nxzt; ++ixzt) {
    int iy = 0;
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int site = ix + Nx * (iy + Ny * izt);
    real_t bc2 = bc[1];
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_yp1(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = bc2 * vt[ivcd];
      }
    }
  }
 }
}

//====================================================================
void mult_domainwall_5din_ym1(real_t *RESTRICT buf_ym,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_by = Nin5bd * (Nst/Ny);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_ym[0:size_by])				\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int idir = 1;
  int Nxzt = Nx * Nz * Nt;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixzt = 0; ixzt < Nxzt; ++ixzt) {
    int iy  = Ny-1;
    int ix  = ixzt % Nx;
    int izt = ixzt / Nx;
    int site = ix + Nx * (iy + Ny * izt);
    real_t ut[NDF];
    load_u(ut, up, site + Nst*idir);
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_ym1(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_zp1(real_t *RESTRICT buf_zp,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_bz = Nin5bd * (Nst/Nz);

#pragma acc data present(up[0:size_u], wp[0:size], buf_zp[0:size_bz])\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {
  int Nxy  = Nx * Ny;
  int Nxyt = Nx * Ny * Nt;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixyt = 0; ixyt < Nxyt; ++ixyt) {
    int iz = 0;
    int ixy = ixyt % Nxy;
    int it  = ixyt / Nxy;
    int site = ixy + Nxy * (iz + Nz * it);
    real_t bc2 = bc[2];
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_zp1(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = bc2 * vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_zm1(real_t *RESTRICT buf_zm,
                              real_t *RESTRICT up, real_t *RESTRICT wp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_bz = Nin5bd * (Nst/Nz);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_zm[0:size_bz])				\
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int idir = 2;
  int Nxy  = Nx * Ny;
  int Nxyt = Nx * Ny * Nt;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixyt = 0; ixyt < Nxyt; ++ixyt) {
    int iz = Nz-1;
    int ixy = ixyt % Nxy;
    int it  = ixyt / Nxy;
    int site = ixy + Nxy * (iz + Nz * it);
    real_t ut[NDF];
    load_u(ut, up, site + Nst*idir);
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_zm1(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_tp1_dirac(
                             real_t *RESTRICT buf_tp,
                             real_t *RESTRICT up, real_t *RESTRICT wp,
                             int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_bt = Nin5bd * (Nst/Nt);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_tp[0:size_bt] ) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {
  int Nxyz = Nx * Ny * Nz;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixyz = 0; ixyz < Nxyz; ++ixyz) {
    int it = 0;
    int site = ixyz + Nxyz * it;
    real_t bc2 = bc[3];
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_tp1_dirac(vt, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = bc2 * vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_tm1_dirac(
                             real_t *RESTRICT buf_tm,
                             real_t *RESTRICT up, real_t *RESTRICT wp,
                             int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * NDIM;
  int size_bt = Nin5bd * (Nst/Nt);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_tm[0:size_bt] ) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

  int idir = 3;
  int Nxyz = Nx * Ny * Nz;

#pragma acc parallel  \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop
  for (int ixyz = 0; ixyz < Nxyz; ++ixyz) {
    int it = Nt-1;
    int site = ixyz + Nxyz * it;
    real_t ut[NDF];
    load_u(ut, up, site + Nst*idir);
    for (int is = 0; is < Ns; ++is) {
      real_t wt[NVCD], vt[NVC * ND2];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }
      mult_wilson_tm1_dirac(vt, ut, wt);
      for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
        buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = vt[ivcd];
      }
    }
  }
 }

}

//====================================================================
void mult_domainwall_5din_xp2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_xp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = Nin5bd * Ny * Nz * Nt;

#pragma acc data present(buf_xp[0:size_bx], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      idir = 0;

      if (ix == Nx-1) {
        real_t ut[NDF];
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

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_xm2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_xm,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bx = Nin5bd * Ny * Nz * Nt;

#pragma acc data present(buf_xm[0:size_bx], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      if (ix == 0) {
        real_t bc2 = bc[0];
        real_t wt[NVC * ND2];
        for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
          wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
        }
        mult_wilson_xm2(vL, wt);

        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	}
      }

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_yp2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_yp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_by = Nin5bd * Nx * Nz * Nt;

#pragma acc data present(buf_yp[0:size_by], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      idir = 1;
      int ixzt = ix + Nx * izt;

      if (iy == Ny-1) {
        real_t ut[NDF];
        load_u(ut, up, site + Nst * idir);
        real_t wt[NVC * ND2];
        for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
          wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
        }
        mult_wilson_yp2(vL, ut, wt);

        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	}
      }

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_ym2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_ym,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_by = Nin5bd * Nx * Nz * Nt;

#pragma acc data present(buf_ym[0:size_by], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int ixzt = ix + Nx * izt;

      if (iy == 0) {
        real_t bc2 = bc[1];
        real_t wt[NVC * ND2];
        for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
          wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
        }
        mult_wilson_ym2(vL, wt);

        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	}
      }

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}

//====================================================================
void mult_domainwall_5din_zp2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_zp,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bz = Nin5bd * Nx * Ny * Nt;

#pragma acc data present(buf_zp[0:size_bz],  \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int idir = 2;
      int ixyt = ixy + Nxy * it;
      if (iz == Nz-1) {
        real_t ut[NDF];
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

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_zm2(real_t *RESTRICT vp, real_t *RESTRICT up,
                              real_t *RESTRICT buf_zm,
                              int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bz = Nin5bd * Nx * Ny * Nt;

#pragma acc data present(buf_zm[0:size_bz], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int ixy  = site % Nxy;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int ixyt = ixy + Nxy * it;

      if (iz == 0) {
        real_t bc2 = bc[2];
        real_t wt[NVC * ND2];
        for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
          wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
        }
        mult_wilson_zm2(vL, wt);

        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	}
      }

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_tp2_dirac(
                             real_t *RESTRICT vp, real_t *RESTRICT up,
                             real_t *RESTRICT buf_tp,
                             int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bt = Nin5bd * Nx * Ny * Nz;

#pragma acc data present(buf_tp[0:size_bt],  \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;
  int Nxyz = Nx * Ny * Nz;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int izt  = site / Nxy;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int idir = 3;
      if (it == Nt-1) {
        real_t ut[NDF];
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

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}
//====================================================================
void mult_domainwall_5din_tm2_dirac(
                             real_t *RESTRICT vp, real_t *RESTRICT up,
                             real_t *RESTRICT buf_tm,
                             int Ns, int *bc, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 4;
  int size_bt = Nin5bd * Nx * Ny * Nz;

#pragma acc data present(buf_tm[0:size_bt], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, bc[0:NDIM])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;
  int Nxyz = Nx * Ny * Nz;

#pragma acc loop
  for (int site = 0; site < Nst; ++site) {
    int izt  = site / Nxy;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      if (it == 0) {
        real_t bc2 = bc[3];
        real_t wt[NVC * ND2];
        for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
          wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
        }
        mult_wilson_tm2_dirac(vL, wt);

        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, (ivcd + NVCD*is), site)] += vL[ivcd];
	}
      }

    }   // is loop
  }   // site loop

 }  // acc parallel
 }  // acc data

}

#endif
//============================================================END=====
