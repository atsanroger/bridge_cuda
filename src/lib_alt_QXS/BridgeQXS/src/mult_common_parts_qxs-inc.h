/*!
      @file    mult_common_parts_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_COMMON_PATRS_QXS_H
#define MULT_COMMON_PATRS_QXS_H


//====================================================================
namespace {
  inline void load_udag_xm(svbool_t pg1, svbool_t pg2,
                           svreal_t& u0, svreal_t& u1,
                           svreal_t& u2, svreal_t& u3,
                           svreal_t& u4, svreal_t& u5,
                           real_t *ux, real_t *un)
  {
    load_vec(pg1, u0, &ux[VLEN * 0 - 1]);
    load_add(pg2, u0, &un[VLEN * 0 + VLENX - 1]);

    load_vec(pg1, u1, &ux[VLEN * 1 - 1]);
    load_add(pg2, u1, &un[VLEN * 1 + VLENX - 1]);

    load_vec(pg1, u2, &ux[VLEN * 2 - 1]);
    load_add(pg2, u2, &un[VLEN * 2 + VLENX - 1]);

    load_vec(pg1, u3, &ux[VLEN * 3 - 1]);
    load_add(pg2, u3, &un[VLEN * 3 + VLENX - 1]);

    load_vec(pg1, u4, &ux[VLEN * 4 - 1]);
    load_add(pg2, u4, &un[VLEN * 4 + VLENX - 1]);

    load_vec(pg1, u5, &ux[VLEN * 5 - 1]);
    load_add(pg2, u5, &un[VLEN * 5 + VLENX - 1]);
  }


  inline void load_udag_ym(svbool_t pg1, svbool_t pg2,
                           svreal_t& u0, svreal_t& u1,
                           svreal_t& u2, svreal_t& u3,
                           svreal_t& u4, svreal_t& u5,
                           real_t *ux, real_t *un)
  {
#if VLENY > 1
    load_vec(pg1, u0, &ux[VLEN * 0 - VLENX]);
    load_add(pg2, u0, &un[VLEN * 0 + VLENX * (VLENY - 1)]);

    load_vec(pg1, u1, &ux[VLEN * 1 - VLENX]);
    load_add(pg2, u1, &un[VLEN * 1 + VLENX * (VLENY - 1)]);

    load_vec(pg1, u2, &ux[VLEN * 2 - VLENX]);
    load_add(pg2, u2, &un[VLEN * 2 + VLENX * (VLENY - 1)]);

    load_vec(pg1, u3, &ux[VLEN * 3 - VLENX]);
    load_add(pg2, u3, &un[VLEN * 3 + VLENX * (VLENY - 1)]);

    load_vec(pg1, u4, &ux[VLEN * 4 - VLENX]);
    load_add(pg2, u4, &un[VLEN * 4 + VLENX * (VLENY - 1)]);

    load_vec(pg1, u5, &ux[VLEN * 5 - VLENX]);
    load_add(pg2, u5, &un[VLEN * 5 + VLENX * (VLENY - 1)]);
#else
    svbool_t pg = set_predicate();
    load_vec(pg, u0, &un[VLEN * 0]);
    load_vec(pg, u1, &un[VLEN * 1]);
    load_vec(pg, u2, &un[VLEN * 2]);
    load_vec(pg, u3, &un[VLEN * 3]);
    load_vec(pg, u4, &un[VLEN * 4]);
    load_vec(pg, u5, &un[VLEN * 5]);
#endif
  }


  inline void load_udag_ym(svreal_t& u0, svreal_t& u1,
                           svreal_t& u2, svreal_t& u3,
                           svreal_t& u4, svreal_t& u5,
                           real_t *ux, real_t *un)
  {
#if VLENY > 1
    shift_vec_yfw(u0, &ux[VLEN * 0], &un[VLEN * 0]);
    shift_vec_yfw(u1, &ux[VLEN * 1], &un[VLEN * 1]);
    shift_vec_yfw(u2, &ux[VLEN * 2], &un[VLEN * 2]);
    shift_vec_yfw(u3, &ux[VLEN * 3], &un[VLEN * 3]);
    shift_vec_yfw(u4, &ux[VLEN * 4], &un[VLEN * 4]);
    shift_vec_yfw(u5, &ux[VLEN * 5], &un[VLEN * 5]);
#else
    svbool_t pg = set_predicate();
    load_vec(pg, u0, &un[VLEN * 0]);
    load_vec(pg, u1, &un[VLEN * 1]);
    load_vec(pg, u2, &un[VLEN * 2]);
    load_vec(pg, u3, &un[VLEN * 3]);
    load_vec(pg, u4, &un[VLEN * 4]);
    load_vec(pg, u5, &un[VLEN * 5]);
#endif
  }


  inline void load_udag_xm2_eo(svbool_t pg1, svbool_t pg3,
                               svreal_t& u0, svreal_t& u1,
                               svreal_t& u2, svreal_t& u3,
                               svreal_t& u4, svreal_t& u5,
                               real_t *ux)
  {
    load_vec(pg3, u0, &ux[VLEN * 0]);
    load_add(pg1, u0, &ux[VLEN * 0 - 1]);

    load_vec(pg3, u1, &ux[VLEN * 1]);
    load_add(pg1, u1, &ux[VLEN * 1 - 1]);

    load_vec(pg3, u2, &ux[VLEN * 2]);
    load_add(pg1, u2, &ux[VLEN * 2 - 1]);

    load_vec(pg3, u3, &ux[VLEN * 3]);
    load_add(pg1, u3, &ux[VLEN * 3 - 1]);

    load_vec(pg3, u4, &ux[VLEN * 4]);
    load_add(pg1, u4, &ux[VLEN * 4 - 1]);

    load_vec(pg3, u5, &ux[VLEN * 5]);
    load_add(pg1, u5, &ux[VLEN * 5 - 1]);
  }


  inline void load_udag_xm_eo(svbool_t pg1, svbool_t pg2, svbool_t pg3,
                              svreal_t& u0, svreal_t& u1,
                              svreal_t& u2, svreal_t& u3,
                              svreal_t& u4, svreal_t& u5,
                              real_t *ux, real_t *un)
  {
    load_vec(pg3, u0, &ux[VLEN * 0]);
    load_add(pg1, u0, &ux[VLEN * 0 - 1]);
    load_add(pg2, u0, &un[VLEN * 0 + VLENX - 1]);

    load_vec(pg3, u1, &ux[VLEN * 1]);
    load_add(pg1, u1, &ux[VLEN * 1 - 1]);
    load_add(pg2, u1, &un[VLEN * 1 + VLENX - 1]);

    load_vec(pg3, u2, &ux[VLEN * 2]);
    load_add(pg1, u2, &ux[VLEN * 2 - 1]);
    load_add(pg2, u2, &un[VLEN * 2 + VLENX - 1]);

    load_vec(pg3, u3, &ux[VLEN * 3]);
    load_add(pg1, u3, &ux[VLEN * 3 - 1]);
    load_add(pg2, u3, &un[VLEN * 3 + VLENX - 1]);

    load_vec(pg3, u4, &ux[VLEN * 4]);
    load_add(pg1, u4, &ux[VLEN * 4 - 1]);
    load_add(pg2, u4, &un[VLEN * 4 + VLENX - 1]);

    load_vec(pg3, u5, &ux[VLEN * 5]);
    load_add(pg1, u5, &ux[VLEN * 5 - 1]);
    load_add(pg2, u5, &un[VLEN * 5 + VLENX - 1]);
  }


  inline void load_udag(svbool_t pg1, svuint_t idx1,
                        svreal_t& u0, svreal_t& u1,
                        svreal_t& u2, svreal_t& u3,
                        svreal_t& u4, svreal_t& u5,
                        real_t *ux, real_t *un)
  {
    shift_vec(pg1, idx1, u0, &ux[VLEN * 0], &un[VLEN * 0]);
    shift_vec(pg1, idx1, u1, &ux[VLEN * 1], &un[VLEN * 1]);
    shift_vec(pg1, idx1, u2, &ux[VLEN * 2], &un[VLEN * 2]);
    shift_vec(pg1, idx1, u3, &ux[VLEN * 3], &un[VLEN * 3]);
    shift_vec(pg1, idx1, u4, &ux[VLEN * 4], &un[VLEN * 4]);
    shift_vec(pg1, idx1, u5, &ux[VLEN * 5], &un[VLEN * 5]);
  }
} // nameless namespace end

#endif
//============================================================END=====
