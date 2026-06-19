/*!
      @file    mult_Doainwall_5din_eo_inv_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2620 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_INV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_INV_ACC_INCLUDED

//#define IDX_MAT(idx1, idx2, Ns)  ((idx1) + (ND2 * Ns) * (idx2)))

//====================================================================
void mult_domainwall_5din_ee_LUinv_dirac_org(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns])
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
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = a * vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
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
        vt[ivcd] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
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
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
      }

    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_LUinv_dirac(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst * NVC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;

    real_t vt[ND], yt[ND], xt[ND];

    int is = 0;
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = e[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
	vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += e[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    vt[0] += a * (xt[0] + xt[2]);
    vt[1] += a * (xt[1] + xt[3]);
    vt[2] += a * (xt[2] + xt[0]);
    vt[3] += a * (xt[3] + xt[1]);

    vt[0] += -0.5 * (yt[0] - yt[2]);
    vt[1] += -0.5 * (yt[1] - yt[3]);
    vt[2] += -0.5 * (yt[2] - yt[0]);
    vt[3] += -0.5 * (yt[3] - yt[1]);

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = a * vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] + vt[2]);
    yt[1] = 0.5 * (vt[1] + vt[3]);
    yt[2] = 0.5 * (vt[2] + vt[0]);
    yt[3] = 0.5 * (vt[3] + vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += - f[is] * yt[id];
      }

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
      }

    }

  } // idx loop

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_LUinv_dirac_alt(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst * NVC * ND2; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivcd1 = idx2_wp % (NVC * ND2);
    int idx_out = idx2_wp / (NVC * ND2);
    int site = idx_in + NWP * idx_out;

    int ivcd2 = ivcd1 + NVC * 2;
    real_t vt1, vt2, yt1, yt2, xt1, xt2;

    int is = 0;
    vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

    yt1 = e[0] * vt1;
    yt2 = e[0] * vt2;

    for (int is = 1; is < Ns-1; ++is) {

      xt1 = vt1;
      xt2 = vt2;

      vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
      vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      vt1 += a * (xt1 + xt2);
      vt2 += a * (xt1 + xt1);

      vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
      vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

      yt1 += e[is] * vt1;
      yt2 += e[is] * vt2;

    }

    is = Ns-1;

    xt1 = vt1;
    xt2 = vt2;

    vt1 = wp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = wp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    vt1 += a * (xt1 + xt2);
    vt2 += a * (xt2 + xt1);

    vt1 += -0.5 * (yt1 - yt2);
    vt2 += -0.5 * (yt2 - yt1);

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];

    vt1 = a * vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
    vt2 = a * vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

    vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)] = vt2;

    yt1 = 0.5 * (vt1 + vt2);
    yt2 = 0.5 * (vt2 + vt1);

    for (int is = Ns-2; is >= 0; --is) {

      xt1 = vt1;
      xt2 = vt2;

      vt1 = vp[IDX2(Nin5, (ivcd1 + NVCD*is), site)];
      vt2 = vp[IDX2(Nin5, (ivcd2 + NVCD*is), site)];

      real_t a = real_t(0.5) * dm[is];

      vt1 += a * (xt1 - xt2);
      vt2 += a * (xt2 - xt1);

      vt1 += - f[is] * yt1;
      vt2 += - f[is] * yt2;

      real_t aa = dpinv[is];

      vt1 *= aa;
      vt2 *= aa;

      vp[IDX2(Nin5, (ivcd1 + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ivcd2 + NVCD * is), site)] = vt2;

    }

  } // idx loop

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_LUdaginv_dirac(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst * NVC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;

    real_t vt[ND], yt[ND], xt[ND];
    //qqqqq
    int is = 0;

    real_t a = dpinv[0];

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = a * wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = f[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is - 1];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += f[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    a = real_t(0.5) * dm[is - 1];

    vt[0] += a * (xt[0] - xt[2]);
    vt[1] += a * (xt[1] - xt[3]);
    vt[2] += a * (xt[2] - xt[0]);
    vt[3] += a * (xt[3] - xt[1]);

    vt[0] += -0.5 * (yt[0] + yt[2]);
    vt[1] += -0.5 * (yt[1] + yt[3]);
    vt[2] += -0.5 * (yt[2] + yt[0]);
    vt[3] += -0.5 * (yt[3] + yt[1]);

    real_t aa = dpinv[is];

    for(int id = 0; id < ND; ++id){
      vt[id] *= aa;
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }
    // Udag_inv completed

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] - vt[2]);
    yt[1] = 0.5 * (vt[1] - vt[3]);
    yt[2] = 0.5 * (vt[2] - vt[0]);
    yt[3] = 0.5 * (vt[3] - vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += -e[is] * yt[id];
      }
      
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD*is), site)] = vt[id];
      }
    }

  } // idx loop end

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_LUdaginv_dirac_org(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;

  int size    = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns])
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
    // Udag_inv completed

    is = Ns-1;

    for(int ivcd = 0; ivcd < NVCD; ++ivcd){
      vt[ivcd] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
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
        vt[ivcd] = vp[IDX2(Nin5, (ivcd + NVCD*is), site)];
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

  } // site loop end

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_inv_dirac_4d(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
			int jd,
                        int Ns, real_t *RESTRICT mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;
  int mat_size  = ND2 * Ns;
  int mat_size2 = mat_size * mat_size;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Ns, Nin5, Nst, mat_inv[0:mat_size2])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst * NC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int idx_out = idx2_wp / NC;
    int site    = idx_in + NWP * idx_out;
   // if(site < Nst){
   //  for (int site = 0; site < Nst; ++site) {

    for (int is1 = 0; is1 < Ns; ++is1) {
      //for (int ic = 0; ic < NC; ++ic) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2, index;
            real_t mat1, mat2;
	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
	    }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
	    }
            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
    }
  }

 }
 }

 }

//====================================================================
void mult_domainwall_5din_ee_inv_dirac_5d(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
			int jd,
                        int Ns, real_t *RESTRICT mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;
  int mat_size  = ND2 * Ns;
  int mat_size2 = mat_size * mat_size;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Ns, Nin5, Nst, mat_inv[0:mat_size2])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst * Ns * NC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site = idx_in + NWP * idx_out;
    // if(site < Nst){
      /*
    for (int idx = 0; idx < Nst_pad * Ns; ++idx) {
      int idx2_wp = idx / NWP;
      int idx_in  = idx % NWP;
      int is1     = idx2_wp % Ns;
      int idx_out = idx2_wp / Ns;
      int site    = idx_in + NWP * idx_out;
      */
      //for (int ic = 0; ic < NC; ++ic) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2, index;
            real_t mat1, mat2;
	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
	    }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
	    }
            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, (0+2*(ic+NC*id1) + NVCD * is1), site)] = vt1r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id1) + NVCD * is1), site)] = vt1i;
          vp[IDX2(Nin5, (0+2*(ic+NC*id2) + NVCD * is1), site)] = vt2r;
          vp[IDX2(Nin5, (1+2*(ic+NC*id2) + NVCD * is1), site)] = vt2i;
	}
      }
    }

 }
 }

}

#endif
//============================================================END=====
