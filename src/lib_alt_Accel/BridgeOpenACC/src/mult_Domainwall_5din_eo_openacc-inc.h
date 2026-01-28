/*!
      @file    mult_Doainwall_5din_eo_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2618 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_ee_5dir_dirac(
       real_t *RESTRICT vp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
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
      }
    }
  }

 }
 }

}

//====================================================================
//void mult_domainwall_5din_eo_5dir_dirac(
void mult_domainwall_5din_eo_5dir_dirac_org(
       real_t *RESTRICT yp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(yp[0:size], wp[0:size]) \
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

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        yp[IDX2(Nin5, (ivcd + NVCD * is), site)]
                 = -0.5 * (b[is] * wt[ivcd] + c[is] * vt[ivcd]);
      }
    }
  }

 }
 }

}

//====================================================================
//void mult_domainwall_5din_eo_5dir_dirac_alt(
void mult_domainwall_5din_eo_5dir_dirac(
       real_t *RESTRICT yp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(yp[0:size], wp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst_pad * NVC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP*idx_out;
    if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t vt1, vt2, vt3, vt4;
      real_t wt1, wt2, wt3, wt4;

      int is_up = (is+1) % Ns;
      real_t Fup   = 0.5;
      if (is == Ns-1) Fup = -0.5 * mq;
      wt1  = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
      wt2  = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
      wt3  = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
      wt4  = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];

      vt1  = Fup * (wt1 - wt3);
      vt2  = Fup * (wt2 - wt4);
      vt3  = Fup * (wt3 - wt1);
      vt4  = Fup * (wt4 - wt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5;
      if (is == 0) Fdn = -0.5 * mq;
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];

      vt1 += Fdn * (wt1+ wt3);
      vt2 += Fdn * (wt2+ wt4);
      vt3 += Fdn * (wt3+ wt1);
      vt4 += Fdn * (wt4+ wt2);

      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];

      yp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)]
                 = -0.5 * (b[is] * wt1 + c[is] * vt1);
      yp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)]
                 = -0.5 * (b[is] * wt2 + c[is] * vt2);
      yp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)]
                 = -0.5 * (b[is] * wt3 + c[is] * vt3);
      yp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)]
                 = -0.5 * (b[is] * wt4 + c[is] * vt4);

    } //is
    } // site < Nst

  } // idx
 }

 }
}

//====================================================================
void mult_domainwall_5din_ee_5dirdag_dirac(
       real_t *RESTRICT vp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD], xt[NVCD];

      real_t B1 = b[is] * (4.0 - M0) + 1.0;
      real_t a1 = -0.5 * b[is];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vt[ivcd] = B1 * wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      int is_up = (is+1) % Ns;
      real_t C1 = c[is_up] * (4.0 - M0) - 1.0;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        xt[ivcd] = C1 * wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
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

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        xt[ivcd] = C2 * wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
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
//void mult_domainwall_5din_eo_5dirdag_dirac(
void mult_domainwall_5din_eo_5dirdag_dirac_org(
       real_t *RESTRICT vp, real_t *RESTRICT yp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], yp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int site = 0; site < Nst; ++site) {

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD], xt[NVCD], yt[NVCD];

      real_t a1 = -0.5 * b[is];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] = a1 * yt[ID3 + ivc];
        vt[ID2 + ivc] = a1 * yt[ID4 + ivc];
        vt[ID3 + ivc] = a1 * yt[ID1 + ivc];
        vt[ID4 + ivc] = a1 * yt[ID2 + ivc];
      }

      int is_up = (is+1) % Ns;
      real_t aup = -0.5 * c[is_up];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        xt[ID1 + ivc] = aup * yt[ID3 + ivc];
        xt[ID2 + ivc] = aup * yt[ID4 + ivc];
        xt[ID3 + ivc] = aup * yt[ID1 + ivc];
        xt[ID4 + ivc] = aup * yt[ID2 + ivc];
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
      real_t adn = -0.5 * c[is_dn];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        yt[ivcd] = yp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        xt[ID1 + ivc] = adn * yt[ID3 + ivc];
        xt[ID2 + ivc] = adn * yt[ID4 + ivc];
        xt[ID3 + ivc] = adn * yt[ID1 + ivc];
        xt[ID4 + ivc] = adn * yt[ID2 + ivc];
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
//void mult_domainwall_5din_eo_5dirdag_dirac_alt(
void mult_domainwall_5din_eo_5dirdag_dirac(
       real_t *RESTRICT vp, real_t *RESTRICT yp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = ceil_nwp(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], yp[0:size]) \
                 copyin(Nst, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst_pad*NVC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP*idx_out;
    if(site<Nst){
      for (int is = 0; is < Ns; ++is) {

      real_t vt1, vt2, vt3, vt4;
      real_t xt1, xt2, xt3, xt4;
      real_t yt1, yt2, yt3, yt4;

      real_t a1 = -0.5 * b[is];

      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)];

      vt1  = a1 * yt3;
      vt2  = a1 * yt4;
      vt3  = a1 * yt1;
      vt4  = a1 * yt2;

      int is_up = (is+1) % Ns;
      real_t aup = -0.5 * c[is_up];

      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_up), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_up), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_up), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_up), site)];

      xt1  = aup * yt3;
      xt2  = aup * yt4;
      xt3  = aup * yt1;
      xt4  = aup * yt2;

      real_t Fup = 0.5;
      if (is == Ns-1) Fup = -0.5 * mq;

      vt1  += Fup * (xt1 + xt3);
      vt2  += Fup * (xt2 + xt4);
      vt3  += Fup * (xt3 + xt1);
      vt4  += Fup * (xt4 + xt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t adn = -0.5 * c[is_dn];

      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD*is_dn), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD*is_dn), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD*is_dn), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD*is_dn), site)];

      xt1  = adn * yt3;
      xt2  = adn * yt4;
      xt3  = adn * yt1;
      xt4  = adn * yt2;

      real_t Fdn   = 0.5;
      if (is == 0) Fdn = -0.5 * mq;

      vt1  += Fdn * (xt1 - xt3);
      vt2  += Fdn * (xt2 - xt4);
      vt3  += Fdn * (xt3 - xt1);
      vt4  += Fdn * (xt4 - xt2);

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      } // is
    } // site < Nst

  } // idx
 }

}
}

//====================================================================
void mult_domainwall_5din_eo_hopb_dirac_4d(
       real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
       int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
       int jgm5)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5 = NVCD * Ns;

  int size = Nin5 * Nst;
  int size_u = NDF * Nst * 2 * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM], ieo, jeo, jgm5)
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

    int keo  = (jeo + iy + iz + it) % 2;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      idir = 0;

      // mult_xp
      {
        int ix2 = (ix + keo) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1 && keo == 1) bc2 = bc[0];
        //if(ix + keo == Nx) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, site + Nst * (ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_xpb(vL, ut, wt);
      }
      // mult_xm
      {
        int ix2 = (ix - 1 + keo + Nx) % Nx;
        int nei  = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0 && keo == 0) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst * (1-ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
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
        load_u(ut, up, site + Nst * (ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_ypb(vL, ut, wt);
      }

      if ((iy > 0) || (do_comm[idir] == 0)) {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc[1];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst * (1-ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
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
        load_u(ut, up, site + Nst * (ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_zpb(vL, ut, wt);
      }

      if ((iz > 0) || (do_comm[idir] == 0)) {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc[2];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst * (1-ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
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
        load_u(ut, up, site + Nst * (ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
        }
        mult_wilson_tpb_dirac(vL, ut, wt);
      }

      if ((it > 0) || (do_comm[idir] == 0)) {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc[3];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst * (1-ieo + 2*idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD*is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), nei)];
          }
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
void mult_domainwall_5din_eo_hopb_dirac_5d(
       real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
       int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
       int jgm5)
{
  int nx  = Nsize[0];
  int ny  = Nsize[1];
  int nz  = Nsize[2];
  int nt  = Nsize[3];
  int nst = nx * ny * nz * nt;
  int Nst_pad = ceil_nwp(nst);

  int Nin5 = NVCD * Ns;

  int size = Nin5 * nst;
  int size_u = NDF * nst * 2 * NDIM;

  real_t *restrict u_up = up;
  real_t *restrict u_dn = up;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size])	\
                 copyin(nst, nx, ny, nz, nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM], ieo, jeo, jgm5)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst_pad * Ns; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;

    int site = idx_in + NWP*idx_out;
    if(site<nst) {

      int nxy  = nx  * ny;
      int nxyz = nxy * nz;

      int ix   = site % nx;
      int iyzt = site / nx;
      int ixy  = site % nxy;
      int iy   = iyzt % ny;
      int izt  = site / nxy;
      int iz   = izt % nz;
      int it   = izt / nz;
      int ixyz = site % nxyz;

      int keo  = (jeo + iy + iz + it) % 2;

      real_t u_0, u_1, u_2, u_3, u_4, u_5;
      real_t u_6, u_7, u_8, u_9, u10, u11;
      real_t u12, u13, u14, u15, u16, u17;
      real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
      real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
      real_t wt1r, wt1i, wt2r, wt2i;

      real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
      real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
      real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
      real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

      real_t *restrict v1 = wp;
      real_t *restrict v2 = vp;
 
#include "inc/mult_Domainwall_eo_openacc2_xyz-inc.h"

#include "inc/mult_Domainwall_eo_openacc2_t_dirac-inc.h"

      v2[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
      v2[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
      v2[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
      v2[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
      v2[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
      v2[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

      v2[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
      v2[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
      v2[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
      v2[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
      v2[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
      v2[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

      v2[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
      v2[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
      v2[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
      v2[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
      v2[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
      v2[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

      v2[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
      v2[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
      v2[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
      v2[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
      v2[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
      v2[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;
     }

   } // site loop
 }  // acc parallel
 }  // acc data

}

//====================================================================
void mult_domainwall_5din_eo_hop1_dirac(
       real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
       real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
       real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
       real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
       real_t *RESTRICT up, real_t *RESTRICT wp,
       int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
       int jgm5)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size   = Nin5 * Nst;
  int size_u = NDF * Nst * 2 * NDIM;

  int size_bx = Nin5bd * ((Ny * Nz * Nt + 1)/2);
  int size_by = Nin5bd * (Nst/Ny);
  int size_bz = Nin5bd * (Nst/Nz);
  int size_bt = Nin5bd * (Nst/Nt);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt] ) \
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:4], ieo, jeo, jgm5)
 {

  if (do_comm[0] > 0) {
    int idir = 0;
    int Nyzt = Ny * Nz * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
#pragma acc loop gang worker vector
    for (int iyzt = 0; iyzt < Nyzt; ++iyzt) {
      int iy = iyzt % Ny;
      int iz = (iyzt/Ny) % Nz;
      int it = iyzt/(Ny * Nz);
      int keo = (jeo + iy + iz + it) % 2;

      if(keo == 1){
        int ix = 0;
        int iyzt2 = iyzt/2;
        int site = ix + Nx * iyzt;
        real_t bc2 = bc[0];
        for (int is = 0; is < Ns; ++is) {
          real_t wt[NVCD], vt[NVC * ND2];
	  if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
            }
	  }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
              wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
              wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
              wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
            }
	  }
          mult_wilson_xp1(vt, wt);
          for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
            buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt2)] = bc2 * vt[ivcd];
          }
	}
      }

      if(keo == 0){
        int iyzt2 = iyzt/2;
        int ix = Nx-1;
        int site = ix + Nx * iyzt;
        real_t ut[NDF];
        load_u(ut, up, site + Nst * (1-ieo + 2*idir));
        for (int is = 0; is < Ns; ++is) {
          real_t wt[NVCD], vt[NVC * ND2];
	  if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
            }
	  }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
              wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
              wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
              wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
            }
	  }
          mult_wilson_xm1(vt, ut, wt);
          for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
            buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt2)] = vt[ivcd];
          }
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
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
        }
        mult_wilson_yp1(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_yp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixzt)] = bc2 * vt[ivcd];
	}
      }

      iy = Ny-1;
      site = ix + Nx * (iy + Ny * izt);
      real_t ut[NDF];
      load_u(ut, up, site + Nst * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
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
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
        }
        mult_wilson_zp1(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_zp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyt)] = bc2 * vt[ivcd];
	}
      }

      iz = Nz-1;
      site = ixy + Nxy * (iz + Nz * it);
      real_t ut[NDF];
      load_u(ut, up, site + Nst * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
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
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
        }
        mult_wilson_tp1_dirac(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_tp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), ixyz)] = bc2 * vt[ivcd];
	}
      }

      it = Nt-1;
      site = ixyz + Nxyz * it;
      real_t ut[NDF];
      load_u(ut, up, site + Nst * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD*is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD*is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD*is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD*is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD*is), site)];
          }
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
void mult_domainwall_5din_eo_hop2_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
      real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
      real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
      real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
      real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
      int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst;
  int size_u  = NDF * Nst * 2 * NDIM;
  int size_bx = Nin5bd * ((Ny * Nz * Nt + 1)/2);
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
                        bc[0:NDIM], do_comm[0:4], ieo, jeo)
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
    int keo  = (jeo + iy + iz + it) % 2;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int opr_any = 0;

      idir = 0;
      if (do_comm[idir] > 0) {

        if(ix == Nx-1 && keo == 1){
          real_t ut[NDF];
          load_u(ut, up, site + Nst * (ieo + 2*idir));
          real_t wt[NVC * ND2];
          int iyzt2 = iyzt/2;
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = buf_xp[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt2)];
          }
          mult_wilson_xp2(vL, ut, wt);
          ++opr_any;
        }

        if(ix == 0 && keo == 0){
          real_t bc2 = bc[0];
          int iyzt2  = iyzt/2;
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVC*ND2; ++ivcd){
            wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, (ivcd + NVC*ND2*is), iyzt2)];
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
          load_u(ut, up, site + Nst * (ieo + 2*idir));
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
          load_u(ut, up, site + Nst * (ieo + 2*idir));
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
          load_u(ut, up, site + Nst * (ieo + 2*idir));
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


#endif
//============================================================END=====
