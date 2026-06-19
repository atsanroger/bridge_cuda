/*!
      @file    mult_Doainwall_5din_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2618 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_5dir_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT yp, real_t *RESTRICT wp,
      real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], yp[0:size], wp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

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

}

//====================================================================
void mult_domainwall_5din_5dirdag_dirac(
       real_t *RESTRICT vp, real_t *RESTRICT yp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], yp[0:size], wp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

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

}


//====================================================================
void mult_domainwall_5din_mult_gm5_dirac(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD], wt[NVCD];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd]  = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] = wt[ID3 + ivc];
        vt[ID2 + ivc] = wt[ID4 + ivc];
        vt[ID3 + ivc] = wt[ID1 + ivc];
        vt[ID4 + ivc] = wt[ID2 + ivc];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
	vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
      }
    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_mult_R(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD];
      int isR = Ns-1 - is;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * isR), site)];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
	vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
      }
    }

  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_mult_gm5R_dirac(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    for (int is = 0; is < Ns; ++is) {

      int isR = Ns-1 - is;
      real_t vt[NVCD], wt[NVCD];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*isR), site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] = wt[ID3 + ivc];
        vt[ID2 + ivc] = wt[ID4 + ivc];
        vt[ID3 + ivc] = wt[ID1 + ivc];
        vt[ID4 + ivc] = wt[ID2 + ivc];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
	vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }

    }
  }  // site loop
 }
 }

}

//====================================================================
void mult_domainwall_5din_clear(real_t *RESTRICT vp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size]) copyin(Nst, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {
    for (int is = 0; is < Ns; ++is) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
	vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = 0.0;
      }
    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_hopb_dirac(
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

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
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

      idir = 1;

      if ((iy < Ny-1) || (do_comm[idir] == 0)) {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc[1];

        real_t ut[NDF];
        load_u(ut, up, site + Nst * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
	}
        mult_wilson_ypb(vL, ut, wt);
      }

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
void mult_domainwall_5din_hop1_dirac(
                  real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                  real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                  real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                  real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                  real_t *RESTRICT up, real_t *RESTRICT wp,
                  int Ns, int *bc, int *Nsize, int *do_comm)
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
  int size_by = Nin5bd * (Nst/Ny);
  int size_bz = Nin5bd * (Nst/Nz);
  int size_bt = Nin5bd * (Nst/Nt);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt] ) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:4])
 {

  if (do_comm[0] > 0) {
    int idir = 0;
    int Nyzt = Ny * Nz * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop gang worker vector
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

      ix = Nx-1;
      site = ix + Nx * iyzt;
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

  } // do_comm[0]

  if (do_comm[1] > 0) {
    int idir = 1;
    int Nxzt = Nx * Nz * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop gang worker vector
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

      iy = Ny-1;
      site = ix + Nx * (iy + Ny * izt);
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

  } // do_comm[1]

  if (do_comm[2] > 0) {
    int idir = 2;
    int Nxy  = Nx * Ny;
    int Nxyt = Nx * Ny * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop gang worker vector
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

      iz = Nz-1;
      site = ixy + Nxy * (iz + Nz * it);
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

  } // do_comm[2]

  if (do_comm[3] > 0) {
    int idir = 3;
    int Nxyz = Nx * Ny * Nz;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop gang worker vector
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

      it = Nt-1;
      site = ixyz + Nxyz * it;
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

  } // do_comm[3]

#pragma acc wait

 } // acc data

  if(do_comm[0] > 0){
    #pragma acc update async host (buf_xp[0:size_bx])
    #pragma acc update async host (buf_xm[0:size_bx])
  }
  if(do_comm[1] > 0){
    #pragma acc update async host (buf_yp[0:size_by])
    #pragma acc update async host (buf_ym[0:size_by])
  }
  if(do_comm[2] > 0){
    #pragma acc update async host (buf_zp[0:size_bz])
    #pragma acc update async host (buf_zm[0:size_bz])
  }
  if(do_comm[3] > 0){
    #pragma acc update async host (buf_tp[0:size_bt])
    #pragma acc update async host (buf_tm[0:size_bt])
  }

 #pragma acc wait

}

//====================================================================
void mult_domainwall_5din_hop2_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
      real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
      real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
      real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
      real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
      int Ns, int *bc, int *Nsize, int *do_comm)
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
  int size_by = Nin5bd * Nx * Nz * Nt;
  int size_bz = Nin5bd * Nx * Ny * Nt;
  int size_bt = Nin5bd * Nx * Ny * Nz;

  if(do_comm[0] > 0){
    #pragma acc update async device (buf_xp[0:size_bx])
    #pragma acc update async device (buf_xm[0:size_bx])
  }
  if(do_comm[1] > 0){
    #pragma acc update async device (buf_yp[0:size_by])
    #pragma acc update async device (buf_ym[0:size_by])
  }
  if(do_comm[2] > 0){
    #pragma acc update async device (buf_zp[0:size_bz])
    #pragma acc update async device (buf_zm[0:size_bz])
  }
  if(do_comm[3] > 0){
    #pragma acc update async device (buf_tp[0:size_bt])
    #pragma acc update async device (buf_tm[0:size_bt])
  }

#pragma acc wait

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:4])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;
  int Nxyz = Nx * Ny * Nz;

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {
    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int opr_any = 0;

      idir = 0;
      if (do_comm[idir] > 0) {

        if (ix == Nx-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst * idir);
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
          }
          mult_wilson_xp2(vL, ut, wt);
          ++opr_any;
        }

        if (ix == 0) {
          real_t bc2 = bc[0];
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt)];
          }
          mult_wilson_xm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 1;
      if (do_comm[idir] > 0) {
        int ixzt = ix + Nx * izt;

        if (iy == Ny-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst * idir);
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
          }
          mult_wilson_yp2(vL, ut, wt);
          ++opr_any;
        }

        if (iy == 0) {
          real_t bc2 = bc[1];
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)];
          }
          mult_wilson_ym2(vL, wt);
          ++opr_any;
        }
      }

      idir = 2;
      if (do_comm[idir] > 0) {
        int ixyt = ixy + Nxy * it;
        if (iz == Nz-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst * idir);
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
          }
          mult_wilson_zp2(vL, ut, wt);
          ++opr_any;
        }

        if (iz == 0) {
          real_t bc2 = bc[2];
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)];
          }
          mult_wilson_zm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 3;
      if (do_comm[idir] > 0) {
        if (it == Nt-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst * idir);
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
          }
          mult_wilson_tp2_dirac(vL, ut, wt);
          ++opr_any;
        }

        if (it == 0) {
          real_t bc2 = bc[3];
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)];
          }
          mult_wilson_tm2_dirac(vL, wt);
          ++opr_any;
        }
      }

      if (opr_any > 0) {
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
void mult_domainwall_5din_L_inv_dirac(
            real_t *RESTRICT vp, real_t *RESTRICT wp,
            int Ns, int *Nsize, real_t *e, real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    real_t vt[NVCD], yt[NVCD], xt[NVCD];

    int is = 0;
    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }
    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      yt[ivcd] = e[0] * vt[ivcd];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        xt[ivcd] = vt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += a * (xt[ID1 + ivc] + xt[ID3 + ivc]);
        vt[ID2 + ivc] += a * (xt[ID2 + ivc] + xt[ID4 + ivc]);
        vt[ID3 + ivc] += a * (xt[ID3 + ivc] + xt[ID1 + ivc]);
        vt[ID4 + ivc] += a * (xt[ID4 + ivc] + xt[ID2 + ivc]);
      }
	
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        yt[ivcd] += e[is] * vt[ivcd];
      }

    }

    is = Ns-1;

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      xt[ivcd] = vt[ivcd];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ID1 + ivc] += a * (xt[ID1 + ivc] + xt[ID3 + ivc]);
      vt[ID2 + ivc] += a * (xt[ID2 + ivc] + xt[ID4 + ivc]);
      vt[ID3 + ivc] += a * (xt[ID3 + ivc] + xt[ID1 + ivc]);
      vt[ID4 + ivc] += a * (xt[ID4 + ivc] + xt[ID2 + ivc]);
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      vt[ID1 + ivc] += -0.5 * (yt[ID1 + ivc] - yt[ID3 + ivc]);
      vt[ID2 + ivc] += -0.5 * (yt[ID2 + ivc] - yt[ID4 + ivc]);
      vt[ID3 + ivc] += -0.5 * (yt[ID3 + ivc] - yt[ID1 + ivc]);
      vt[ID4 + ivc] += -0.5 * (yt[ID4 + ivc] - yt[ID2 + ivc]);
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
    }

  } // site loop

 }
 }

}


//====================================================================
void mult_domainwall_5din_U_inv_dirac(
             real_t *RESTRICT vp, real_t *RESTRICT wp,
             int Ns, int *Nsize, real_t *f, real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, f[0:Ns-1], dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    real_t vt[NVCD], yt[NVCD], xt[NVCD];

    int is = Ns-1;

    real_t a = dpinv[Ns-1];

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = a * wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      yt[ID1 + ivc] = 0.5 * (vt[ID1 + ivc] + vt[ID3 + ivc]);
      yt[ID2 + ivc] = 0.5 * (vt[ID2 + ivc] + vt[ID4 + ivc]);
      yt[ID3 + ivc] = 0.5 * (vt[ID3 + ivc] + vt[ID1 + ivc]);
      yt[ID4 + ivc] = 0.5 * (vt[ID4 + ivc] + vt[ID2 + ivc]);
    }

    for (int is = Ns-2; is >= 0; --is) {

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        xt[ivcd] = vt[ivcd];
      }
	
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += a * (xt[ID1 + ivc] - xt[ID3 + ivc]);
        vt[ID2 + ivc] += a * (xt[ID2 + ivc] - xt[ID4 + ivc]);
        vt[ID3 + ivc] += a * (xt[ID3 + ivc] - xt[ID1 + ivc]);
        vt[ID4 + ivc] += a * (xt[ID4 + ivc] - xt[ID2 + ivc]);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] += - f[is] * yt[ivcd];
      }

      real_t aa = dpinv[is];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] *= aa;
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }

    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_Ldag_inv_dirac(
            real_t *RESTRICT vp, real_t *RESTRICT wp,
            int Ns, int *Nsize, real_t *e, real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    real_t vt[NVCD], yt[NVCD], xt[NVCD];

    int is = Ns-1;
    
    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
    }

    for(int ivc = 0; ivc < NVC; ++ivc){
      yt[ID1 + ivc] = 0.5 * (vt[ID1 + ivc] - vt[ID3 + ivc]);
      yt[ID2 + ivc] = 0.5 * (vt[ID2 + ivc] - vt[ID4 + ivc]);
      yt[ID3 + ivc] = 0.5 * (vt[ID3 + ivc] - vt[ID1 + ivc]);
      yt[ID4 + ivc] = 0.5 * (vt[ID4 + ivc] - vt[ID2 + ivc]);
    }

    for (int is = Ns-2; is >= 0; --is) {

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        xt[ivcd] = vt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += a * (xt[ID1 + ivc] + xt[ID3 + ivc]);
        vt[ID2 + ivc] += a * (xt[ID2 + ivc] + xt[ID4 + ivc]);
        vt[ID3 + ivc] += a * (xt[ID3 + ivc] + xt[ID1 + ivc]);
        vt[ID4 + ivc] += a * (xt[ID4 + ivc] + xt[ID2 + ivc]);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] += -e[is] * yt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }
    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_Udag_inv_dirac(
            real_t *RESTRICT vp, real_t *RESTRICT wp,
            int Ns, int *Nsize, real_t *f, real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, f[0:Ns-1], dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    real_t vt[NVCD], yt[NVCD], xt[NVCD];

    int is = 0;

    real_t a = dpinv[0];

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = a * wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      yt[ivcd] = f[0] * vt[ivcd];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        xt[ivcd] = vt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is - 1];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += a * (xt[ID1 + ivc] - xt[ID3 + ivc]);
        vt[ID2 + ivc] += a * (xt[ID2 + ivc] - xt[ID4 + ivc]);
        vt[ID3 + ivc] += a * (xt[ID3 + ivc] - xt[ID1 + ivc]);
        vt[ID4 + ivc] += a * (xt[ID4 + ivc] - xt[ID2 + ivc]);
      }
	
      real_t aa = dpinv[is];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] *= aa;
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        yt[ivcd] += f[is] * vt[ivcd];
      }

    }

    is = Ns-1;
    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      xt[ivcd] = vt[ivcd];
    }

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    a = real_t(0.5) * dm[is - 1];

      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += a * (xt[ID1 + ivc] - xt[ID3 + ivc]);
        vt[ID2 + ivc] += a * (xt[ID2 + ivc] - xt[ID4 + ivc]);
        vt[ID3 + ivc] += a * (xt[ID3 + ivc] - xt[ID1 + ivc]);
        vt[ID4 + ivc] += a * (xt[ID4 + ivc] - xt[ID2 + ivc]);
      }
      for(int ivc = 0; ivc < NVC; ++ivc){
        vt[ID1 + ivc] += -0.5 * (yt[ID1 + ivc] + yt[ID3 + ivc]);
        vt[ID2 + ivc] += -0.5 * (yt[ID2 + ivc] + yt[ID4 + ivc]);
        vt[ID3 + ivc] += -0.5 * (yt[ID3 + ivc] + yt[ID1 + ivc]);
        vt[ID4 + ivc] += -0.5 * (yt[ID4 + ivc] + yt[ID2 + ivc]);
      }

      real_t aa = dpinv[is];
      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vt[ivcd] *= aa;
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[ivcd];
      }

  } // site loop

 }
 }

}


#endif
//============================================================END=====
