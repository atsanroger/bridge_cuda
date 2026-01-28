// afopr_Wilson_mult-inc.h

namespace {

//====================================================================
inline void check_setup()
{
  if(VLEN2 < 2){
    vout.crucial("VLEN2 = %d is too small for this implementation\n",
                 VLEN2);
    exit(EXIT_FAILURE);
  }

}
//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xp1(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  REALTYPE vt[NVCD];

  load_vec1(vt, v1, 0, NCD);
  set_sp2_xp1(buf, vt, Nc);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xp2(REALTYPE *v2, REALTYPE *u,
                            REALTYPE *buf, int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  set_th(vt1, REALTYPE(0.0), VLEN*NC);
  set_th(vt2, REALTYPE(0.0), VLEN*NC);

  for(int ic = 0; ic < NC; ++ic){
    vt1[VLEN-2 + VLEN*ic] = buf[2*ic  ];
    vt1[VLEN-1 + VLEN*ic] = buf[2*ic+1];
    vt2[VLEN-2 + VLEN*ic] = buf[2*ic   + NVC];
    vt2[VLEN-1 + VLEN*ic] = buf[2*ic+1 + NVC];
  }

  load_vec(ut, u, NDF2);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_xp(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xpb(REALTYPE *v2, REALTYPE *u,
                            REALTYPE *v1, int Nc)
{
  REALTYPE vt[VLEN*NCD], wt[VLEN*NCD];
  REALTYPE ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(wt, v1, NCD);
  set_sp2_xp(vt1, vt2, wt, Nc);

  load_vec(ut, u,  NDF2);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_xp(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xm1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN*NC], wt2[VLEN*NC];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_xm(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(&wt1[ic2], &ut[ic2], vt1, Nc);
    mult_udagv(&wt2[ic2], &ut[ic2], vt2, Nc);
  }

  for(int ic = 0; ic < NC; ++ic){
    buf[2*ic  ]       = wt1[VLEN-2 + VLEN*ic];
    buf[2*ic+1]       = wt1[VLEN-1 + VLEN*ic];
    buf[2*ic   + NVC] = wt2[VLEN-2 + VLEN*ic];
    buf[2*ic+1 + NVC] = wt2[VLEN-1 + VLEN*ic];
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xmb(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], wt[VLEN*NCD];
  REALTYPE ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN*NC], wt2[VLEN*NC];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_xm(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(&wt1[ic2], &ut[ic2], vt1, Nc);
    mult_udagv(&wt2[ic2], &ut[ic2], vt2, Nc);
  }

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    set_sp4_xm(&v2[ic2], &wt1[ic2], &wt2[ic2], Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xm2(REALTYPE *v2, REALTYPE *buf, int Nc)
{
  REALTYPE yt1[NVC], yt2[NVC];

  for(int ivc = 0; ivc < NVC; ++ivc){
    yt1[ivc] = buf[ivc];
    yt2[ivc] = buf[ivc + NVC];
  }

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = 2 * ic;
    set_sp4_xm1(&v2[VLEN*ic], 0, &yt1[ic2], &yt2[ic2], Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xm2(REALTYPE *v2, REALTYPE *buf,
                            int site, int ibf, int Nc)
{
  REALTYPE yt1[NVC], yt2[NVC];

  int iv = VLEN * NCD * site;
  int ib = NVC * ND2 * ibf;

  int ib2 = ib + NVC;
  for(int ivc = 0; ivc < NVC; ++ivc){
    yt1[ivc] = buf[ivc + ib];
    yt2[ivc] = buf[ivc + ib2];
  }

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = 2 * ic;
    int icv = VLEN * ic + iv;
    set_sp4_xm1(&v2[icv], 0, &yt1[ic2], &yt2[ic2], Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_yp1(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  REALTYPE vt[VLEN*NCD];

  load_vec(vt, v1, NCD);
  set_sp2_yp(&buf[0], &buf[VLEN*NC], vt, Nc);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_yp2(REALTYPE *v2, REALTYPE *u, REALTYPE *buf,
                            int Nc)
{
  REALTYPE ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(ut, u, NDF2);

  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN*NC], NC);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_yp(&v2[VLEN * ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ypb(REALTYPE *v2, REALTYPE *u,
                            REALTYPE *v1, int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_yp(vt1, vt2, vt, Nc);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_yp(&v2[VLEN * ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ym1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);
  set_sp2_ym(vt1, vt2, v1, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(&buf[ic2],           &ut[ic2], vt1, Nc);
    mult_udagv(&buf[ic2 + VLEN*NC], &ut[ic2], vt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ym2(REALTYPE *v2, REALTYPE *buf, int Nc)
{
  REALTYPE wt1[VLEN], wt2[VLEN];

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    load_vec(wt1, &buf[ic2],           1);
    load_vec(wt2, &buf[ic2 + VLEN*NC], 1);
    set_sp4_ym(&v2[VLEN * ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ymb(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);
  set_sp2_ym(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(wt1, &ut[ic2], vt1, Nc);
    mult_udagv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_ym(&v2[ic2], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zp1(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  REALTYPE vt[VLEN*NCD];

  load_vec(vt, v1, NCD);
  set_sp2_zp(&buf[0], &buf[VLEN * NC], vt, Nc);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zp2(REALTYPE *v2, REALTYPE *u, REALTYPE *buf,
                            int Nc)
{
  REALTYPE ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(ut, u, NDF2);
  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN * NC], NC);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_zp(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zpb(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_zp(vt1, vt2, vt, Nc);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_zp(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zm1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);
  set_sp2_zm(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(&buf[ic2],         &ut[ic2], vt1, Nc);
    mult_udagv(&buf[ic2+VLEN*NC], &ut[ic2], vt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zm2(REALTYPE *v2, REALTYPE *buf, int Nc)
{
  REALTYPE wt1[VLEN], wt2[VLEN];

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    load_vec(wt1, &buf[ic2],  1);
    load_vec(wt2, &buf[ic2 + VLEN * NC], 1);
    set_sp4_zm(&v2[ic2], wt1, wt2, Nc);
  }
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zmb(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);
  set_sp2_zm(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(wt1, &ut[ic2], vt1, Nc);
    mult_udagv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_zm(&v2[ic2], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tp1_dirac(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  REALTYPE vt[VLEN*NCD];

  load_vec(vt, v1, NCD);
  set_sp2_tp_dirac(&buf[0], &buf[VLEN*NC], vt, Nc);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tp2_dirac(REALTYPE *v2, REALTYPE *u, REALTYPE *buf,
                                  int Nc)
{
  REALTYPE ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(ut, u, NDF2);

  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN*NC], NC);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_tp_dirac(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tpb_dirac(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_tp_dirac(vt1, vt2, vt, Nc);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_tp_dirac(&v2[VLEN*ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tm1_dirac(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);
  set_sp2_tm_dirac(vt1, vt2, vt, Nc);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(&buf[ic2],           &ut[ic2], vt1, Nc);
    mult_udagv(&buf[ic2 + VLEN*NC], &ut[ic2], vt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tm2_dirac(REALTYPE *v2, REALTYPE *buf, int Nc)
{
  REALTYPE wt1[VLEN], wt2[VLEN];

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    load_vec(wt1, &buf[ic2],           1);
    load_vec(wt2, &buf[ic2 + VLEN*NC], 1);
    set_sp4_tm_dirac(&v2[ic2], wt1, wt2, Nc);
  }
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tmb_dirac(REALTYPE *v2, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  REALTYPE vt[VLEN*NCD], ut[VLEN*NDF2];
  REALTYPE vt1[VLEN*NC], vt2[VLEN*NC];
  REALTYPE wt1[VLEN], wt2[VLEN];

  load_vec(vt, v1, NCD);
  load_vec(ut, u,  NDF2);

  set_sp2_tm_dirac(vt1, vt2, vt, Nc);

  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(wt1, &ut[ic2], vt1, Nc);
    mult_udagv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_tm_dirac(&v2[ic2], wt1, wt2, Nc);
  }

}

//====================================================================

} // nameless namespace end
