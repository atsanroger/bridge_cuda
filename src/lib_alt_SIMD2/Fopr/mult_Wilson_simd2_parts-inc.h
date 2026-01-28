/*!
      @file    mult_Wilson_simd2_parts-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef MULT_WILSON_SIMD2_PARTS_H
#define MULT_WILSON_SIMD2_PARTS_H

namespace {

//====================================================================
inline void check_setup()
{
  /*
  if(VLEN2 < 2){
    vout.crucial("VLEN2 = %d is too small for this implementation\n",
                 VLEN2);
    exit(EXIT_FAILURE);
  }
  */
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
inline void mult_wilson_xp2(Vsimd_t *v2,
                             REALTYPE *u, REALTYPE *buf, int Nc)
{
  Vsimd_t ut[NDF2];
  Vsimd_t vt1[NC], vt2[NC];
  Vsimd_t wt1, wt2;
  REALTYPE zero[VLEN];

  set_th(zero, REALTYPE(0.0), VLEN);

  for(int ic = 0; ic < NC; ++ic){
    shift_vec1_bw(&vt1[ic], zero, &buf[2*ic], 1);
    shift_vec1_bw(&vt2[ic], zero, &buf[2*ic + NVC], 1);
  }

  load_vec(ut, u, NDF2);
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_xp(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xpb(Vsimd_t *v2,
                            REALTYPE *u, REALTYPE *v1, int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_xp(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_xp(&v2[ic], wt1, wt2, Nc);
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
  /*
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_xm(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    // need masked store. the following code does not work.
    save_vec(&buf[VLEN*ic], &wt1, 1);
    save_vec(&buf[VLEN*(ic+NC)], &wt2, 1);
  }
  */

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xm2(Vsimd_t *v2, REALTYPE *buf, int Nc)
{
  REALTYPE zero[VLEN];
  set_th(zero, REALTYPE(0.0), VLEN);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    shift_vec1_fw(&wt1, zero, &buf[2*ic], 1);
    shift_vec1_fw(&wt2, zero, &buf[2*ic + NVC], 1);
    set_sp4_xm(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_xmb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_xm(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    set_sp4_xm(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_yp1(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_yp(vt1, vt2, v1, Nc);
  save_vec(&buf[0], vt1, NC);
  save_vec(&buf[VLEN*NC], vt2, NC);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_yp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN*NC], NC);

  Vsimd_t ut[NDF2];
  load_vec(ut, u, NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_yp(&v2[ic], wt1, wt2, Nc);
  }

 }

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ypb(Vsimd_t *v2, REALTYPE *u,
                            REALTYPE *v1, int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_yp(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u, NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_yp(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ym1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_ym(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    save_vec(&buf[VLEN*ic], &wt1, 1);
    save_vec(&buf[VLEN*(ic+NC)], &wt2, 1);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ym2(Vsimd_t *v2, REALTYPE *buf, int Nc)
{
  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    load_vec(&wt1, &buf[VLEN*ic],        1);
    load_vec(&wt2, &buf[VLEN*(ic + NC)], 1);
    set_sp4_ym(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_ymb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_ym(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    set_sp4_ym(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zp1(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_zp(vt1, vt2, v1, Nc);
  save_vec(&buf[0], vt1, NC);
  save_vec(&buf[VLEN*NC], vt2, NC);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zp2(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN * NC], NC);

  Vsimd_t ut[NDF2];
  load_vec(ut, u, NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_zp(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zpb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_zp(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_zp(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zm1(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_zm(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    save_vec(&buf[VLEN*ic], &wt1, 1);
    save_vec(&buf[VLEN*(ic+NC)], &wt2, 1);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zm2(Vsimd_t *v2, REALTYPE *buf, int Nc)
{
  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    load_vec(&wt1, &buf[VLEN*ic],        1);
    load_vec(&wt2, &buf[VLEN*(ic + NC)], 1);
    set_sp4_zm(&v2[ic], wt1, wt2, Nc);
  }
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_zmb(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                            int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_zm(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    set_sp4_zm(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tp1_dirac(REALTYPE *buf, REALTYPE *v1, int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_tp_dirac(vt1, vt2, v1, Nc);
  save_vec(&buf[0], vt1, NC);
  save_vec(&buf[VLEN*NC], vt2, NC);
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tp2_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *buf,
                                  int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_vec(vt1, &buf[0],  NC);
  load_vec(vt2, &buf[VLEN*NC], NC);

  Vsimd_t ut[NDF2];
  load_vec(ut, u, NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_tp_dirac(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tpb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_tp_dirac(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = NC * ic;
    mult_uv(wt1, &ut[ic2], vt1, Nc);
    mult_uv(wt2, &ut[ic2], vt2, Nc);
    set_sp4_tp_dirac(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tm1_dirac(REALTYPE *buf, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_tm_dirac(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    save_vec(&buf[VLEN*ic], &wt1, 1);
    save_vec(&buf[VLEN*(ic+NC)], &wt2, 1);
  }

}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tm2_dirac(Vsimd_t *v2, REALTYPE *buf, int Nc)
{
  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    int ic2 = VLEN * ic;
    load_vec(&wt1, &buf[ic2],           1);
    load_vec(&wt2, &buf[ic2 + VLEN*NC], 1);
    set_sp4_tm_dirac(&v2[ic], wt1, wt2, Nc);
  }
}

//====================================================================
template <typename REALTYPE>
inline void mult_wilson_tmb_dirac(Vsimd_t *v2, REALTYPE *u, REALTYPE *v1,
                                  int Nc)
{
  Vsimd_t vt1[NC], vt2[NC];
  load_set_sp2_tm_dirac(vt1, vt2, v1, Nc);

  Vsimd_t ut[NDF2];
  load_vec(ut, u,  NDF2);

  Vsimd_t wt1, wt2;
  for(int ic = 0; ic < NC; ++ic){
    mult_udagv(wt1, &ut[ic], vt1, Nc);
    mult_udagv(wt2, &ut[ic], vt2, Nc);
    set_sp4_tm_dirac(&v2[ic], wt1, wt2, Nc);
  }

}

//====================================================================

} // nameless namespace end

#endif
