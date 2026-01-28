/*!
      @file    mult_Clover_coarse_parts_qxs-inc.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (Mon, 20 Mar 2023) $
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_CLOVER_COARSE_PARTS_QXS_H
#define MULT_CLOVER_COARSE_PARTS_QXS_H

namespace {
//====================================================================
  inline void accum_mult_u_i(svreal_t& out_e1r, svreal_t& out_e1i,
                             svreal_t& out_e2r, svreal_t& out_e2i,
                             const real_t *__restrict__ in,
                             const real_t *__restrict__ u0,
                             int i, const int ncol)
  {
    const int nh = ncol / 2;
    enum velem
    {
      e1r, e1i, e2r, e2i
    };
    svbool_t pg = set_predicate();

    const real_t *u = u0 + VLEN * 4 * i * ncol;
    for (int j = 0; j < nh; j++) {
      svreal_t     u_e1r, u_e1i, u_e2r, u_e2i;
      svreal_t     in_e1r, in_e1i, in_e2r, in_e2i;
      const real_t *uij = u + VLEN * 4 * j;
      const real_t *inj = in + VLEN * 4 * j;
      load_vec(pg, u_e1r, uij);
      load_vec(pg, u_e1i, uij + VLEN);
      load_vec(pg, u_e2r, uij + 2 * VLEN);
      load_vec(pg, u_e2i, uij + 3 * VLEN);
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      //  out1 += (u[2*j  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      ymax_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      axpy_vec(pg, out_e1i, u_e1i, in_e1r);

      // out1 += u[2*j+1) * in[2*j+1];   // ch:-+
      axpy_vec(pg, out_e1r, u_e2r, in_e2r);
      ymax_vec(pg, out_e1r, u_e2i, in_e2i);
      axpy_vec(pg, out_e1i, u_e2r, in_e2i);
      axpy_vec(pg, out_e1i, u_e2i, in_e2r);

      svreal_t     uu_e1r, uu_e1i, uu_e2r, uu_e2i;
      const real_t *uuij = u + VLEN * 4 * (j + nh);
      load_vec(pg, uu_e1r, uuij);
      load_vec(pg, uu_e1i, uuij + VLEN);
      load_vec(pg, uu_e2r, uuij + 2 * VLEN);
      load_vec(pg, uu_e2i, uuij + 3 * VLEN);

      //  out2 += (u[2*j  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e2r, uu_e1r, in_e1r);
      ymax_vec(pg, out_e2r, uu_e1i, in_e1i);
      axpy_vec(pg, out_e2i, uu_e1r, in_e1i);
      axpy_vec(pg, out_e2i, uu_e1i, in_e1r);

      // out2 += u[2*j+1) * in[2*j+1];   // ch:-+
      axpy_vec(pg, out_e2r, uu_e2r, in_e2r);
      ymax_vec(pg, out_e2r, uu_e2i, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2r, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2i, in_e2r);
    } // j
  }


//====================================================================
  inline void accum_mult_u_i(real_t *__restrict__ out,
                             const real_t *__restrict__ in,
                             const real_t *__restrict__ u0,
                             int i, const int ncol, const svbool_t pgin)
  {
    const int nh = ncol / 2;
    svbool_t  pg = set_predicate();
    svreal_t  out_e1r, out_e1i, out_e2r, out_e2i;
    load_vec(pgin, out_e1r, out);
    load_vec(pgin, out_e1i, out + VLEN);
    load_vec(pgin, out_e2r, out + 2 * VLEN);
    load_vec(pgin, out_e2i, out + 3 * VLEN);

    const real_t *u = u0 + VLEN * 4 * i * ncol;
    for (int j = 0; j < nh; j++) {
      svreal_t     u_e1r, u_e1i, u_e2r, u_e2i;
      svreal_t     in_e1r, in_e1i, in_e2r, in_e2i;
      const real_t *uij = u + VLEN * 4 * j;
      const real_t *inj = in + VLEN * 4 * j;
      load_vec(pg, u_e1r, uij);
      load_vec(pg, u_e1i, uij + VLEN);
      load_vec(pg, u_e2r, uij + 2 * VLEN);
      load_vec(pg, u_e2i, uij + 3 * VLEN);
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      //  out1 += (u[2*j  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      ymax_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      axpy_vec(pg, out_e1i, u_e1i, in_e1r);

      // out1 += u[2*j+1) * in[2*j+1];   // ch:-+
      axpy_vec(pg, out_e1r, u_e2r, in_e2r);
      ymax_vec(pg, out_e1r, u_e2i, in_e2i);
      axpy_vec(pg, out_e1i, u_e2r, in_e2i);
      axpy_vec(pg, out_e1i, u_e2i, in_e2r);

      svreal_t     uu_e1r, uu_e1i, uu_e2r, uu_e2i;
      const real_t *uuij = u + VLEN * 4 * (j + nh);
      load_vec(pg, uu_e1r, uuij);
      load_vec(pg, uu_e1i, uuij + VLEN);
      load_vec(pg, uu_e2r, uuij + 2 * VLEN);
      load_vec(pg, uu_e2i, uuij + 3 * VLEN);

      //  out2 += (u[2*j  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e2r, uu_e1r, in_e1r);
      ymax_vec(pg, out_e2r, uu_e1i, in_e1i);
      axpy_vec(pg, out_e2i, uu_e1r, in_e1i);
      axpy_vec(pg, out_e2i, uu_e1i, in_e1r);

      // out2 += u[2*j+1) * in[2*j+1];   // ch:-+
      axpy_vec(pg, out_e2r, uu_e2r, in_e2r);
      ymax_vec(pg, out_e2r, uu_e2i, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2r, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2i, in_e2r);
    } // j
    save_vec(pg, out, out_e1r);
    save_vec(pg, out + VLEN, out_e1i);
    save_vec(pg, out + 2 * VLEN, out_e2r);
    save_vec(pg, out + 3 * VLEN, out_e2i);
  }


//====================================================================
  inline void accum_mult_u(real_t *out, const real_t *in, const real_t *u0,
                           const int ncol)
  { // simple implementation
    const int nh = ncol / 2;

    for (int i = 0; i < nh; i++) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_u_i(outi, in, u0, i, ncol, set_predicate());
    }
  }


//====================================================================
  inline void accum_mult_udag_i(svbool_t pg,
                                svreal_t& out_e1r, svreal_t& out_e1i,
                                svreal_t& out_e2r, svreal_t& out_e2i,
                                const real_t *__restrict__ in,
                                const real_t *__restrict__ u0,
                                const int i, const int ncol)
  {
    const int dof2 = ncol * ncol;
    const int nh   = ncol / 2;

    for (int j = 0; j < nh; j++) {
      const real_t *u = u0 + VLEN * 4 * j * ncol;
      svreal_t     u_e1r, u_e1i, u_e2r, u_e2i;
      const real_t *uji = u + VLEN * 4 * i;
      load_vec(pg, u_e1r, uji);
      load_vec(pg, u_e1i, uji + VLEN);
      load_vec(pg, u_e2r, uji + 2 * VLEN);
      load_vec(pg, u_e2i, uji + 3 * VLEN);

      svreal_t     in_e1r, in_e1i, in_e2r, in_e2i;
      const real_t *inj = in + VLEN * 4 * j;
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      axpy_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      ymax_vec(pg, out_e1i, u_e1i, in_e1r);

      ymax_vec(pg, out_e2r, u_e2r, in_e1r);
      ymax_vec(pg, out_e2r, u_e2i, in_e1i);
      ymax_vec(pg, out_e2i, u_e2r, in_e1i);
      axpy_vec(pg, out_e2i, u_e2i, in_e1r);

      svreal_t     uu_e1r, uu_e1i, uu_e2r, uu_e2i;
      const real_t *uuji = u + VLEN * (2 * ncol + 4 * i);
      load_vec(pg, uu_e1r, uuji);
      load_vec(pg, uu_e1i, uuji + VLEN);
      load_vec(pg, uu_e2r, uuji + 2 * VLEN);
      load_vec(pg, uu_e2i, uuji + 3 * VLEN);

      ymax_vec(pg, out_e1r, uu_e1r, in_e2r);
      ymax_vec(pg, out_e1r, uu_e1i, in_e2i);
      ymax_vec(pg, out_e1i, uu_e1r, in_e2i);
      axpy_vec(pg, out_e1i, uu_e1i, in_e2r);

      axpy_vec(pg, out_e2r, uu_e2r, in_e2r);
      axpy_vec(pg, out_e2r, uu_e2i, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2r, in_e2i);
      ymax_vec(pg, out_e2i, uu_e2i, in_e2r);
    }
  }


//====================================================================
  inline void accum_mult_udag_i(real_t *__restrict__ out,
                                const real_t *__restrict__ in,
                                const real_t *__restrict__ u0,
                                int i, const int ncol,
                                const svbool_t pgin)
  {
    const int dof2 = ncol * ncol;
    const int nh   = ncol / 2;

    svbool_t pg = set_predicate();
    svreal_t out_e1r, out_e1i, out_e2r, out_e2i;
    load_vec(pgin, out_e1r, out);
    load_vec(pgin, out_e1i, out + VLEN);
    load_vec(pgin, out_e2r, out + 2 * VLEN);
    load_vec(pgin, out_e2i, out + 3 * VLEN);

    for (int j = 0; j < nh; j++) {
      const real_t *u = u0 + VLEN * 4 * j * ncol;

      svreal_t     u_e1r, u_e1i, u_e2r, u_e2i;
      const real_t *uji = u + VLEN * 4 * i;
      load_vec(pg, u_e1r, uji);
      load_vec(pg, u_e1i, uji + VLEN);
      load_vec(pg, u_e2r, uji + 2 * VLEN);
      load_vec(pg, u_e2i, uji + 3 * VLEN);

      svreal_t     in_e1r, in_e1i, in_e2r, in_e2i;
      const real_t *inj = in + VLEN * 4 * j;
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      axpy_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      ymax_vec(pg, out_e1i, u_e1i, in_e1r);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      ymax_vec(pg, out_e2r, u_e2r, in_e1r);
      ymax_vec(pg, out_e2r, u_e2i, in_e1i);
      ymax_vec(pg, out_e2i, u_e2r, in_e1i);
      axpy_vec(pg, out_e2i, u_e2i, in_e1r);

      svreal_t     uu_e1r, uu_e1i, uu_e2r, uu_e2i;
      const real_t *uuji = u + VLEN * (2 * ncol + 4 * i);
      load_vec(pg, uu_e1r, uuji);
      load_vec(pg, uu_e1i, uuji + VLEN);
      load_vec(pg, uu_e2r, uuji + 2 * VLEN);
      load_vec(pg, uu_e2i, uuji + 3 * VLEN);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      ymax_vec(pg, out_e1r, uu_e1r, in_e2r);
      ymax_vec(pg, out_e1r, uu_e1i, in_e2i);
      ymax_vec(pg, out_e1i, uu_e1r, in_e2i);
      axpy_vec(pg, out_e1i, uu_e1i, in_e2r);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      axpy_vec(pg, out_e2r, uu_e2r, in_e2r);
      axpy_vec(pg, out_e2r, uu_e2i, in_e2i);
      axpy_vec(pg, out_e2i, uu_e2r, in_e2i);
      ymax_vec(pg, out_e2i, uu_e2i, in_e2r);
    }

    save_vec(pg, out, out_e1r);
    save_vec(pg, out + VLEN, out_e1i);
    save_vec(pg, out + 2 * VLEN, out_e2r);
    save_vec(pg, out + 3 * VLEN, out_e2i);
  }


//====================================================================
  inline void accum_mult_udag_i(real_t *out, const real_t *in, const real_t *u0,
                                const int i, const int ncol)
  {
    accum_mult_udag_i(out, in, u0, i, ncol, set_predicate());
  }


//====================================================================
  inline void accum_mult_udag_xm_i(svbool_t& pg1, svbool_t& pg2,
                                   real_t *__restrict__ out,
                                   real_t *__restrict__ in,
                                   real_t *u1, real_t *u2,
                                   int i, const int ncol)
  {
    //  u1 and u2 may overlap
    const int nh = ncol / 2;
    svbool_t  pg = set_predicate();
    svreal_t  out_e1r, out_e1i, out_e2r, out_e2i;
    load_vec(pg, out_e1r, out);
    load_vec(pg, out_e1i, out + VLEN);
    load_vec(pg, out_e2r, out + 2 * VLEN);
    load_vec(pg, out_e2i, out + 3 * VLEN);

    for (int j = 0; j < nh; j++) {
      real_t *uu1 = u1 + VLEN * 4 * j * ncol;
      real_t *uu2 = u2 + VLEN * 4 * j * ncol;

      // shifted gauge field: for ch -*
      svreal_t u_e1r, u_e1i, u_e2r, u_e2i;
      shift_vec_xfw(pg1, pg2, u_e1r, &uu1[VLEN * (4 * i)],
                    &uu2[VLEN * (4 * i)]);
      shift_vec_xfw(pg1, pg2, u_e1i, &uu1[VLEN * (4 * i + 1)],
                    &uu2[VLEN * (4 * i + 1)]);
      shift_vec_xfw(pg1, pg2, u_e2r, &uu1[VLEN * (4 * i + 2)],
                    &uu2[VLEN * (4 * i + 2)]);
      shift_vec_xfw(pg1, pg2, u_e2i, &uu1[VLEN * (4 * i + 3)],
                    &uu2[VLEN * (4 * i + 3)]);

      const real_t *inj = in + VLEN * 4 * j;

      svreal_t in_e1r, in_e1i, in_e2r, in_e2i;
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      axpy_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      ymax_vec(pg, out_e1i, u_e1i, in_e1r);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      ymax_vec(pg, out_e2r, u_e2r, in_e1r);
      ymax_vec(pg, out_e2r, u_e2i, in_e1i);
      ymax_vec(pg, out_e2i, u_e2r, in_e1i);
      axpy_vec(pg, out_e2i, u_e2i, in_e1r);

      // shifted gauge field: for ch +*
      shift_vec_xfw(pg1, pg2, u_e1r, &uu1[VLEN * (2 * ncol + 4 * i)],
                    &uu2[VLEN * (2 * ncol + 4 * i)]);
      shift_vec_xfw(pg1, pg2, u_e1i, &uu1[VLEN * (2 * ncol + 4 * i + 1)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 1)]);
      shift_vec_xfw(pg1, pg2, u_e2r, &uu1[VLEN * (2 * ncol + 4 * i + 2)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 2)]);
      shift_vec_xfw(pg1, pg2, u_e2i, &uu1[VLEN * (2 * ncol + 4 * i + 3)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 3)]);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      ymax_vec(pg, out_e1r, u_e1r, in_e2r);
      ymax_vec(pg, out_e1r, u_e1i, in_e2i);
      ymax_vec(pg, out_e1i, u_e1r, in_e2i);
      axpy_vec(pg, out_e1i, u_e1i, in_e2r);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      axpy_vec(pg, out_e2r, u_e2r, in_e2r);
      axpy_vec(pg, out_e2r, u_e2i, in_e2i);
      axpy_vec(pg, out_e2i, u_e2r, in_e2i);
      ymax_vec(pg, out_e2i, u_e2i, in_e2r);
    }

    save_vec(pg, out, out_e1r);
    save_vec(pg, out + VLEN, out_e1i);
    save_vec(pg, out + 2 * VLEN, out_e2r);
    save_vec(pg, out + 3 * VLEN, out_e2i);
  }


//====================================================================
  inline void accum_mult_udag_ym_i(svbool_t& pg1, svbool_t& pg2,
                                   real_t *__restrict__ out,
                                   real_t *__restrict__ in,
                                   real_t *u1, real_t *u2,
                                   int i, const int ncol)
  {
    const int nh = ncol / 2;
    svbool_t  pg = set_predicate();
    svreal_t  out_e1r, out_e1i, out_e2r, out_e2i;
    load_vec(pg, out_e1r, out);
    load_vec(pg, out_e1i, out + VLEN);
    load_vec(pg, out_e2r, out + 2 * VLEN);
    load_vec(pg, out_e2i, out + 3 * VLEN);

    for (int j = 0; j < nh; j++) {
      real_t *uu1 = u1 + VLEN * 4 * j * ncol;
      real_t *uu2 = u2 + VLEN * 4 * j * ncol;

      // shifted gauge field: for ch -*
      svreal_t u_e1r, u_e1i, u_e2r, u_e2i;
      shift_vec_yfw(pg1, pg2, u_e1r, &uu1[VLEN * (4 * i)],
                    &uu2[VLEN * (4 * i)]);
      shift_vec_yfw(pg1, pg2, u_e1i, &uu1[VLEN * (4 * i + 1)],
                    &uu2[VLEN * (4 * i + 1)]);
      shift_vec_yfw(pg1, pg2, u_e2r, &uu1[VLEN * (4 * i + 2)],
                    &uu2[VLEN * (4 * i + 2)]);
      shift_vec_yfw(pg1, pg2, u_e2i, &uu1[VLEN * (4 * i + 3)],
                    &uu2[VLEN * (4 * i + 3)]);

      real_t *inj = in + VLEN * 4 * j;

      svreal_t in_e1r, in_e1i, in_e2r, in_e2i;
      load_vec(pg, in_e1r, inj);
      load_vec(pg, in_e1i, inj + VLEN);
      load_vec(pg, in_e2r, inj + 2 * VLEN);
      load_vec(pg, in_e2i, inj + 3 * VLEN);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      axpy_vec(pg, out_e1r, u_e1r, in_e1r);
      axpy_vec(pg, out_e1r, u_e1i, in_e1i);
      axpy_vec(pg, out_e1i, u_e1r, in_e1i);
      ymax_vec(pg, out_e1i, u_e1i, in_e1r);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      ymax_vec(pg, out_e2r, u_e2r, in_e1r);
      ymax_vec(pg, out_e2r, u_e2i, in_e1i);
      ymax_vec(pg, out_e2i, u_e2r, in_e1i);
      axpy_vec(pg, out_e2i, u_e2i, in_e1r);

      // shifted gauge field: for ch +*
      shift_vec_yfw(pg1, pg2, u_e1r, &uu1[VLEN * (2 * ncol + 4 * i)],
                    &uu2[VLEN * (2 * ncol + 4 * i)]);
      shift_vec_yfw(pg1, pg2, u_e1i, &uu1[VLEN * (2 * ncol + 4 * i + 1)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 1)]);
      shift_vec_yfw(pg1, pg2, u_e2r, &uu1[VLEN * (2 * ncol + 4 * i + 2)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 2)]);
      shift_vec_yfw(pg1, pg2, u_e2i, &uu1[VLEN * (2 * ncol + 4 * i + 3)],
                    &uu2[VLEN * (2 * ncol + 4 * i + 3)]);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      ymax_vec(pg, out_e1r, u_e1r, in_e2r);
      ymax_vec(pg, out_e1r, u_e1i, in_e2i);
      ymax_vec(pg, out_e1i, u_e1r, in_e2i);
      axpy_vec(pg, out_e1i, u_e1i, in_e2r);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      axpy_vec(pg, out_e2r, u_e2r, in_e2r);
      axpy_vec(pg, out_e2r, u_e2i, in_e2i);
      axpy_vec(pg, out_e2i, u_e2r, in_e2i);
      ymax_vec(pg, out_e2i, u_e2i, in_e2r);
    }

    save_vec(pg, out, out_e1r);
    save_vec(pg, out + VLEN, out_e1i);
    save_vec(pg, out + 2 * VLEN, out_e2r);
    save_vec(pg, out + 3 * VLEN, out_e2i);
  }


//====================================================================
  inline void accum_mult_udag(real_t *out, const real_t *in, const real_t *u0,
                              const int ncol)
  {
    // c.c. of U(X,i,s;Y,j,t)  = s t U(Y,j,t;X,i,s)
    // X,Y: coordiante
    // i,j: label for testvectors
    // s,t: +-1 for chirality
    // mult of
    //  sum_{j,t} U(X,i,s; X-mu,j,t) f(X-mu,j,t)
    //  = sum_{j,t}  U(X-mu,j,t; X,i,s)^* s t f(X-mu,j,t)
    //  [^* is complex conjugate]

    const int nh = ncol / 2;
    for (int i = 0; i < nh; ++i) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_udag_i(outi, in, u0, i, ncol);
    }
  }


//====================================================================
  inline void set_mult_u(real_t *out,
                         const real_t *in, const real_t *u0,
                         const int ncol)
  {
    const int nh = ncol / 2;
    //  svbool_t pgfalse = svpfalse();
    svbool_t pgfalse = set_predicate_false();

    for (int i = 0; i < nh; i++) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_u_i(outi, in, u0, i, ncol, pgfalse);
    }
  }


//====================================================================
  inline void set_mult_udag(real_t *out,
                            const real_t *in, const real_t *u,
                            const int ncol)
  {
    const int nh = ncol / 2;
    //  svbool_t pgfalse = svpfalse_b();
    svbool_t pgfalse = set_predicate_false();

    for (int i = 0; i < nh; i++) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_udag_i(outi, in, u, i, ncol, pgfalse);
    }
  }


//====================================================================
  inline void accum_buf(real_t *__restrict__ out,
                        real_t *__restrict__ in, const int ncol)
  {
    svbool_t pg = set_predicate();
    for (int i = 0; i < 2 * ncol; i++) { // 2 for complex
      svreal_t vin, vout;
      load_vec(pg, vin, &in[VLEN * i]);
      load_vec(pg, vout, &out[VLEN * i]);
      add_vec(pg, vout, vin);
      save_vec(pg, &out[VLEN * i], vout);
    }
  }


//====================================================================
  inline void copy_buf(real_t *__restrict__ out,
                       real_t *__restrict__ in, const int ncol)
  {
    svbool_t pg = set_predicate();

    for (int i = 0; i < VLEN * 2 * ncol; i += VLEN) { // 2 for complex
      svreal_t vinout;
      load_vec(pg, vinout, &in[i]);
      save_vec(pg, &out[i], vinout);
    }
  }


//====================================================================
  inline void mult_coarse_xp1(svbool_t& pg2, svint_t& svidx,
                              real_t *__restrict__ buf, real_t *__restrict__ in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg2, vt1, &in[VLEN * (4 * i)]);
      load_vec(pg2, vt2, &in[VLEN * (4 * i + 1)]);
      load_vec(pg2, vt3, &in[VLEN * (4 * i + 2)]);
      load_vec(pg2, vt4, &in[VLEN * (4 * i + 3)]);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i)], vt1, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 1)], vt2, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 2)], vt3, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 3)], vt4, svidx);
    }
  }


//====================================================================
  inline void mult_coarse_xp2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              real_t *__restrict__ out,
                              real_t *__restrict__ u0,
                              real_t *__restrict__ in0,
                              real_t *__restrict__ buf,
                              const int ncol, real_t *__restrict__ work)
  {
    svbool_t  pg = set_predicate();
    const int nh = ncol / 2;

    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg1, vt1, &in0[VLEN * (4 * i) + 1]);
      load_vec(pg1, vt2, &in0[VLEN * (4 * i + 1) + 1]);
      load_vec(pg1, vt3, &in0[VLEN * (4 * i + 2) + 1]);
      load_vec(pg1, vt4, &in0[VLEN * (4 * i + 3) + 1]);
      load_add_gather(pg2, vt1, &buf[VLENY * (4 * i)], svidx);
      load_add_gather(pg2, vt2, &buf[VLENY * (4 * i + 1)], svidx);
      load_add_gather(pg2, vt3, &buf[VLENY * (4 * i + 2)], svidx);
      load_add_gather(pg2, vt4, &buf[VLENY * (4 * i + 3)], svidx);
      save_vec(pg, &work[VLEN * (4 * i)], vt1);
      save_vec(pg, &work[VLEN * (4 * i + 1)], vt2);
      save_vec(pg, &work[VLEN * (4 * i + 2)], vt3);
      save_vec(pg, &work[VLEN * (4 * i + 3)], vt4);
    }

    accum_mult_u(out, work, u0, ncol);
  }


//====================================================================
  inline void mult_coarse_xpb(svbool_t& pg1, svbool_t& pg2,
                              real_t *__restrict__ out,
                              real_t *__restrict__ u,
                              real_t *in1, real_t *in2,
                              const int ncol, real_t *__restrict__ work)
  {
    real_t *in = work;

    int      ncol2 = 2 * ncol;
    svbool_t pg    = set_predicate();
    for (int i = 0; i < ncol2; ++i) {
      svreal_t wt;
      shift_vec_xbw(pg1, pg2, wt, &in1[VLEN * i], &in2[VLEN * i]);
      save_vec(pg, &in[VLEN * i], wt);
    }
    accum_mult_u(out, in, u, ncol);
  }


//====================================================================
  inline void mult_coarse_xm1(svbool_t& pg2, svint_t& svidx,
                              real_t *buf, real_t *u, real_t *in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      clear_vec(pg2, vt1);
      clear_vec(pg2, vt2);
      clear_vec(pg2, vt3);
      clear_vec(pg2, vt4);
      accum_mult_udag_i(pg2, vt1, vt2, vt3, vt4, in, u, i, ncol);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i)], vt1, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 1)], vt2, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 2)], vt3, svidx);
      save_vec_scatter(pg2, &buf[VLENY * (4 * i + 3)], vt4, svidx);
    }
  }


//====================================================================
  inline void mult_coarse_xm2(svbool_t& pg1, svbool_t& pg2, svint_t& svidx,
                              real_t *out, real_t *u0, real_t *in0,
                              real_t *buf, const int ncol)
  {
    svbool_t  pg = set_predicate();
    const int nh = ncol / 2;

    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg, vt1, &out[VLEN * (4 * i)]);
      load_vec(pg, vt2, &out[VLEN * (4 * i + 1)]);
      load_vec(pg, vt3, &out[VLEN * (4 * i + 2)]);
      load_vec(pg, vt4, &out[VLEN * (4 * i + 3)]);
      accum_mult_udag_i(pg1, vt1, vt2, vt3, vt4,
                        &in0[-1], &u0[-1], i, ncol);

      svreal_t wt1, wt2, wt3, wt4;
      load_vec_gather(pg2, wt1, &buf[VLENY * (4 * i)], svidx);
      load_vec_gather(pg2, wt2, &buf[VLENY * (4 * i + 1)], svidx);
      load_vec_gather(pg2, wt3, &buf[VLENY * (4 * i + 2)], svidx);
      load_vec_gather(pg2, wt4, &buf[VLENY * (4 * i + 3)], svidx);
      add_vec(pg2, vt1, wt1);
      add_vec(pg2, vt2, wt2);
      add_vec(pg2, vt3, wt3);
      add_vec(pg2, vt4, wt4);

      save_vec(pg, &out[VLEN * (4 * i)], vt1);
      save_vec(pg, &out[VLEN * (4 * i + 1)], vt2);
      save_vec(pg, &out[VLEN * (4 * i + 2)], vt3);
      save_vec(pg, &out[VLEN * (4 * i + 3)], vt4);
    }
  }


//====================================================================
  inline void mult_coarse_xmb(svbool_t& pg1, svbool_t& pg2,
                              real_t *__restrict__ out,
                              real_t *u1, real_t *u2,
                              real_t *in1, real_t *in2,
                              const int ncol, real_t *__restrict__ work)
  {
    real_t *in = work;

    int      ncol2 = 2 * ncol;
    svbool_t pg    = set_predicate();
    for (int i = 0; i < ncol2; ++i) {
      svreal_t wt;
      shift_vec_xfw(pg1, pg2, wt, &in1[VLEN * i], &in2[VLEN * i]);
      save_vec(pg, &in[VLEN * i], wt);
    }

    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_udag_xm_i(pg1, pg2, outi, in, u1, u2, i, ncol);
    }
  }


//====================================================================
  inline void mult_coarse_yp1(svbool_t& pg2,
                              real_t *__restrict__ buf,
                              real_t *__restrict__ in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg2, vt1, &in[VLEN * (4 * i)]);
      load_vec(pg2, vt2, &in[VLEN * (4 * i + 1)]);
      load_vec(pg2, vt3, &in[VLEN * (4 * i + 2)]);
      load_vec(pg2, vt4, &in[VLEN * (4 * i + 3)]);
      save_vec(pg2, &buf[VLENX * (4 * i)], vt1);
      save_vec(pg2, &buf[VLENX * (4 * i + 1)], vt2);
      save_vec(pg2, &buf[VLENX * (4 * i + 2)], vt3);
      save_vec(pg2, &buf[VLENX * (4 * i + 3)], vt4);
    }
  }


//====================================================================
  inline void mult_coarse_yp2(svbool_t& pg1, svbool_t& pg2,
                              real_t *__restrict__ out,
                              real_t *__restrict__ u0,
                              real_t *__restrict__ in0,
                              real_t *__restrict__ buf,
                              const int ncol, real_t *work)
  {
    svbool_t  pg = set_predicate();
    const int nh = ncol / 2;

    int offset = -VLENX * (VLENY - 1);
    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg1, vt1, &in0[VLEN * (4 * i) + VLENX]);
      load_vec(pg1, vt2, &in0[VLEN * (4 * i + 1) + VLENX]);
      load_vec(pg1, vt3, &in0[VLEN * (4 * i + 2) + VLENX]);
      load_vec(pg1, vt4, &in0[VLEN * (4 * i + 3) + VLENX]);
      load_add(pg2, vt1, &buf[offset + VLENY * (4 * i)]);
      load_add(pg2, vt2, &buf[offset + VLENY * (4 * i + 1)]);
      load_add(pg2, vt3, &buf[offset + VLENY * (4 * i + 2)]);
      load_add(pg2, vt4, &buf[offset + VLENY * (4 * i + 3)]);
      save_vec(pg, &work[VLEN * (4 * i)], vt1);
      save_vec(pg, &work[VLEN * (4 * i + 1)], vt2);
      save_vec(pg, &work[VLEN * (4 * i + 2)], vt3);
      save_vec(pg, &work[VLEN * (4 * i + 3)], vt4);
    }

    accum_mult_u(out, work, u0, ncol);
  }


//====================================================================
  inline void mult_coarse_ypb(svbool_t& pg1, svbool_t& pg2,
                              real_t *__restrict__ out,
                              real_t *__restrict__ u,
                              real_t *in1, real_t *in2,
                              const int ncol, real_t *work)
  {
    real_t *in = work;

    int      ncol2 = 2 * ncol;
    svbool_t pg    = set_predicate();
    for (int i = 0; i < ncol2; ++i) {
      svreal_t wt;
      shift_vec_ybw(pg1, pg2, wt, &in1[VLEN * i], &in2[VLEN * i]);
      save_vec(pg, &in[VLEN * i], wt);
    }
    const int nh = ncol / 2;
    for (int i = 0; i < nh; ++i) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_u_i(outi, in, u, i, ncol, set_predicate());
    }
  }


//====================================================================
  inline void mult_coarse_ym1(svbool_t& pg2,
                              real_t *__restrict__ buf,
                              real_t *__restrict__ u,
                              real_t *__restrict__ in,
                              const int ncol)
  {
    const int nh = ncol / 2;

    int offset = -VLENX * (VLENY - 1);
    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      clear_vec(pg2, vt1);
      clear_vec(pg2, vt2);
      clear_vec(pg2, vt3);
      clear_vec(pg2, vt4);
      accum_mult_udag_i(pg2, vt1, vt2, vt3, vt4, in, u, i, ncol);
      save_vec(pg2, &buf[offset + VLENX * (4 * i)], vt1);
      save_vec(pg2, &buf[offset + VLENX * (4 * i + 1)], vt2);
      save_vec(pg2, &buf[offset + VLENX * (4 * i + 2)], vt3);
      save_vec(pg2, &buf[offset + VLENX * (4 * i + 3)], vt4);
    }
  }


//====================================================================
  inline void mult_coarse_ym2(svbool_t& pg1, svbool_t& pg2,
                              real_t *out, real_t *u0, real_t *in0,
                              real_t *buf, const int ncol)
  {
    svbool_t  pg = set_predicate();
    const int nh = ncol / 2;

    for (int i = 0; i < nh; ++i) {
      svreal_t vt1, vt2, vt3, vt4;
      load_vec(pg, vt1, &out[VLEN * (4 * i)]);
      load_vec(pg, vt2, &out[VLEN * (4 * i + 1)]);
      load_vec(pg, vt3, &out[VLEN * (4 * i + 2)]);
      load_vec(pg, vt4, &out[VLEN * (4 * i + 3)]);
      accum_mult_udag_i(pg1, vt1, vt2, vt3, vt4,
                        &in0[-VLENX], &u0[-VLENX], i, ncol);

      svreal_t wt1, wt2, wt3, wt4;
      load_vec(pg2, wt1, &buf[VLENX * (4 * i)]);
      load_vec(pg2, wt2, &buf[VLENX * (4 * i + 1)]);
      load_vec(pg2, wt3, &buf[VLENX * (4 * i + 2)]);
      load_vec(pg2, wt4, &buf[VLENX * (4 * i + 3)]);
      add_vec(pg2, vt1, wt1);
      add_vec(pg2, vt2, wt2);
      add_vec(pg2, vt3, wt3);
      add_vec(pg2, vt4, wt4);

      save_vec(pg, &out[VLEN * (4 * i)], vt1);
      save_vec(pg, &out[VLEN * (4 * i + 1)], vt2);
      save_vec(pg, &out[VLEN * (4 * i + 2)], vt3);
      save_vec(pg, &out[VLEN * (4 * i + 3)], vt4);
    }
  }


//====================================================================
  inline void mult_coarse_ymb(svbool_t& pg1, svbool_t& pg2,
                              real_t *__restrict__ out,
                              real_t *u1, real_t *u2,
                              real_t *in1, real_t *in2,
                              const int ncol, real_t *__restrict__ work)
  {
    real_t *in = work;

    int      ncol2 = 2 * ncol;
    svbool_t pg    = set_predicate();
    for (int i = 0; i < ncol2; ++i) {
      svreal_t wt;
      shift_vec_yfw(pg1, pg2, wt, &in1[VLEN * i], &in2[VLEN * i]);
      save_vec(pg, &in[VLEN * i], wt);
    }

    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      real_t *outi = out + VLEN * 4 * i;
      accum_mult_udag_ym_i(pg1, pg2, outi, in, u1, u2, i, ncol);
    }
  }


//====================================================================
  inline void mult_coarse_zp1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


//====================================================================
  inline void mult_coarse_zp2(real_t *out, real_t *u, real_t *buf,
                              const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


//====================================================================
  inline void mult_coarse_zpb(real_t *out, real_t *u, real_t *in,
                              const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
  }


//====================================================================
  inline void mult_coarse_zm1(real_t *out, real_t *u, real_t *buf,
                              const int ncol)
  {
    set_mult_udag(out, buf, u, ncol);
  }


//====================================================================
  inline void mult_coarse_zm2(real_t *out, real_t *buf, const int ncol)
  {
    accum_buf(out, buf, ncol);
  }


//====================================================================
  inline void mult_coarse_zmb(real_t *out, real_t *u, real_t *in,
                              const int ncol)
  {
    accum_mult_udag(out, in, u, ncol);
  }


//====================================================================
  inline void mult_coarse_tp1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


//====================================================================
  inline void mult_coarse_tp2(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


//====================================================================
  inline void mult_coarse_tpb(real_t *out, real_t *u, real_t *in,
                              const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
  }


//====================================================================
  inline void mult_coarse_tm1(real_t *out, real_t *u, real_t *buf,
                              const int ncol)
  {
    set_mult_udag(out, buf, u, ncol);
  }


//====================================================================
  inline void mult_coarse_tm2(real_t *out, real_t *buf, const int ncol)
  {
    accum_buf(out, buf, ncol);
  }


//====================================================================
  inline void mult_coarse_tmb(real_t *out, real_t *u, real_t *in,
                              const int ncol)
  {
    accum_mult_udag(out, in, u, ncol);
  }


//====================================================================
} // nameless namespace end

#endif
