/*!
      @file    afopr_Clover_coarse-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"

#define AFOPR_CLOVER_COARSE_TIMER

#ifdef AFOPR_CLOVER_COARSE_TIMER
#include "lib/Tools/timer.h"
#define TIMER_mult_start               timer_mult->start();
#define TIMER_mult_stop                timer_mult->stop();
#define TIMER_pack_start               timer_pack->start();
#define TIMER_pack_stop                timer_pack->stop();
#define TIMER_bulk_start               timer_bulk->start();
#define TIMER_bulk_stop                timer_bulk->stop();
#define TIMER_boundary_start           timer_boundary->start();
#define TIMER_boundary_stop            timer_boundary->stop();
#define TIMER_comm_start               timer_comm->start();
#define TIMER_comm_stop                timer_comm->stop();
#define TIMER_comm_recv_wait_start     timer_comm_recv_wait->start();
#define TIMER_comm_recv_wait_stop      timer_comm_recv_wait->stop();
#define TIMER_comm_send_wait_start     timer_comm_send_wait->start();
#define TIMER_comm_send_wait_stop      timer_comm_send_wait->stop();
#define TIMER_comm_recv_start_start    timer_comm_recv_start->start();
#define TIMER_comm_recv_start_stop     timer_comm_recv_start->stop();
#define TIMER_comm_send_start_start    timer_comm_send_start->start();
#define TIMER_comm_send_start_stop     timer_comm_send_start->stop();
#define TIMER_comm_test_all_start      timer_comm_test_all->start();
#define TIMER_comm_test_all_stop       timer_comm_test_all->stop();
#define TIMER_clear_start              timer_clear->start();
#define TIMER_clear_stop               timer_clear->stop();
#else
#define TIMER_mult_start
#define TIMER_mult_stop
#define TIMER_pack_start
#define TIMER_pack_stop
#define TIMER_bulk_start
#define TIMER_bulk_stop
#define TIMER_boundary_start
#define TIMER_boundary_stop
#define TIMER_comm_start
#define TIMER_comm_stop
#define TIMER_comm_recv_wait_start
#define TIMER_comm_recv_wait_stop
#define TIMER_comm_send_wait_start
#define TIMER_comm_send_wait_stop
#define TIMER_comm_recv_start_start
#define TIMER_comm_recv_start_stop
#define TIMER_comm_send_start_start
#define TIMER_comm_send_start_stop
#define TIMER_comm_test_all_start
#define TIMER_comm_test_all_stop
#define TIMER_clear_start
#define TIMER_clear_stop
#endif


//====================================================================
namespace {
#ifndef QXS_DATA_ALIGNMENT
  constexpr int alignment = 256;
#else
  constexpr int alignment = QXS_DATA_ALIGNMENT;
#endif

  //====================================================================
  static inline void accum_mult_u_i(Vsimd_t *out,
                                    real_t *in,
                                    real_t *u0,
                                    int i,
                                    const int ncol)
  {
    const int nh = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };

    real_t *u = u0 + VLEN * 4 * i * ncol;
    for (int j = 0; j < nh; j++) {
      Vsimd_t vu[4], vin[4];
      load_vec(vu, &u[VLEN * 4 * j], 4);
      load_vec(vin, &in[VLEN * 4 * j], 4);

      //  out1 += (u[2*j  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out1 += u[2*j+1) * in[2*j+1];   // ch:-+
      add_dot_vec(&out[e1r], &vu[e2r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e2r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e2i], &vin[e2r], 1);

      load_vec(vu, &u[VLEN * 4 * (j + nh)], 4);
      //  out2 += (u[2*j  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e2r], &vu[e1r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e1r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e1i], &vin[e1r], 1);

      // out2 += u[2*j+1) * in[2*j+1];   // ch:-+
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    } // j
  }


  //====================================================================
  static inline void accum_mult_u_yp_i(Vsimd_t *out,
                                       real_t *in1, real_t *in2,
                                       real_t *u0,
                                       int i,
                                       const int ncol)
  {
    const int nh = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };

    real_t *u = u0 + VLEN * 4 * ncol * i;
    for (int j = 0; j < nh; j++) {
      Vsimd_t vu[4], vin[4];
      load_vec(vu, &u[VLEN * 4 * j], 4);
      // shifted input vector
      real_t in[VLEN * 4];
      shift_vec2_ybw(in, &in1[VLEN * 4 * j], &in2[VLEN * 4 * j], 4);


      //  out1 += (u[2*j  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out1 += u[2*j+1) * in[2*j+1];   // ch:-+
      add_dot_vec(&out[e1r], &vu[e2r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e2r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e2i], &vin[e2r], 1);

      load_vec(vu, &u[VLEN * 4 * (j + nh)], 4);
      //  out2 += (u[2*j  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e2r], &vu[e1r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e1r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e1i], &vin[e1r], 1);

      // out2 += u[2*j+1) * in[2*j+1];   // ch:-+
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    } // j
  }


  //====================================================================
  static inline void accum_mult_u(real_t *out,
                                  real_t *in,
                                  real_t *u0,
                                  const int ncol)
  { // simple implementation
    /*
    for(int i=0; i<ncol; i++){
      for(int j=0; j<ncol; j++){
        out[i]+=u0[i*ncol+j]*in[j];
      }
    }
    */
    const int nh = ncol / 2;

    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_u_i(tmp, in, u0, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  //====================================================================
  static inline void accum_mult_udag_i(Vsimd_t *out,
                                       real_t *in,
                                       real_t *u0,
                                       int i,
                                       const int ncol)
  {
    const int dof2 = ncol * ncol;
    const int nh   = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };
    for (int j = 0; j < nh; j++) {
      //        const complex_t *u = u0 + 2 * j * ncol;
      real_t *u = u0 + VLEN * 4 * j * ncol;

      Vsimd_t vu[4], vin[4];
      load_vec(vu, &u[VLEN * 4 * i], 4);
      load_vec(vin, &in[VLEN * 4 * j], 4);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      add_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      sub_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      sub_dot_vec(&out[e2r], &vu[e2r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e1i], 1);
      sub_dot_vec(&out[e2i], &vu[e2r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e1r], 1);

      load_vec(vu, &u[VLEN * (2 * ncol + 4 * i)], 4);
      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      sub_dot_vec(&out[e1r], &vu[e1r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e2i], 1);
      sub_dot_vec(&out[e1i], &vu[e1r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e2r], 1);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      add_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      sub_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    }
  }


  //====================================================================
  static inline void accum_mult_udag_xm_i(Vsimd_t *out,
                                          real_t *in,
                                          real_t *u1, real_t *u2,
                                          int i,
                                          const int ncol)
  {
    const int nh = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };
    for (int j = 0; j < nh; j++) {
      //        const complex_t *u = u0 + 2 * j * ncol;
      real_t *uu1 = u1 + VLEN * 4 * j * ncol;
      real_t *uu2 = u2 + VLEN * 4 * j * ncol;

      // work area for field shifting
      real_t u[VLEN * 4];

      // vector variables: treat shifted varaibles
      Vsimd_t vu[4], vin[4];

      // shifted input vector
      load_vec(vin, &in[VLEN * 4 * j], 4);

      // shifted gauge field: for ch -*
      shift_vec2_xfw(u, &uu1[VLEN * 4 * i], &uu2[VLEN * 4 * i], 4);
      load_vec(vu, u, 4);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      add_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      sub_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      sub_dot_vec(&out[e2r], &vu[e2r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e1i], 1);
      sub_dot_vec(&out[e2i], &vu[e2r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e1r], 1);

      // shifted gauge field: for ch +*
      shift_vec2_xfw(u, &uu1[VLEN * (2 * ncol + 4 * i)], &uu2[VLEN * (2 * ncol + 4 * i)], 4);
      load_vec(vu, u, 4);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      sub_dot_vec(&out[e1r], &vu[e1r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e2i], 1);
      sub_dot_vec(&out[e1i], &vu[e1r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e2r], 1);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      add_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      sub_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    }
  }


  //====================================================================
  static inline void accum_mult_udag_ym_i(Vsimd_t *out,
                                          real_t *in1, real_t *in2,
                                          real_t *u1, real_t *u2,
                                          int i,
                                          const int ncol)
  {
    const int nh = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };
    for (int j = 0; j < nh; j++) {
      //        const complex_t *u = u0 + 2 * j * ncol;
      real_t *uu1 = u1 + VLEN * 4 * j * ncol;
      real_t *uu2 = u2 + VLEN * 4 * j * ncol;

      // work area for field shifting
      real_t u[VLEN * 4];
      real_t in[VLEN * 4];

      // vector variables: treat shifted varaibles
      Vsimd_t vu[4], vin[4];

      // shifted input vector
      shift_vec2_yfw(in, &in1[VLEN * 4 * j], &in2[VLEN * 4 * j], 4);
      load_vec(vin, in, 4);

      // shifted gauge field: for ch -*
      shift_vec2_yfw(u, &uu1[VLEN * 4 * i], &uu2[VLEN * 4 * i], 4);
      load_vec(vu, u, 4);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      add_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      sub_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      sub_dot_vec(&out[e2r], &vu[e2r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e1i], 1);
      sub_dot_vec(&out[e2i], &vu[e2r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e1r], 1);

      // shifted gauge field: for ch +*
      shift_vec2_yfw(u, &uu1[VLEN * (2 * ncol + 4 * i)], &uu2[VLEN * (2 * ncol + 4 * i)], 4);
      load_vec(vu, u, 4);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      sub_dot_vec(&out[e1r], &vu[e1r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e2i], 1);
      sub_dot_vec(&out[e1i], &vu[e1r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e2r], 1);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      add_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      sub_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    }
  }


  //====================================================================
  static inline void accum_mult_udag_ym_i(Vsimd_t *out,
                                          real_t *in,
                                          real_t *u1, real_t *u2,
                                          int i,
                                          const int ncol)
  {
    const int nh = ncol / 2;

    enum velem
    {
      e1r, e1i, e2r, e2i
    };
    for (int j = 0; j < nh; j++) {
      //        const complex_t *u = u0 + 2 * j * ncol;
      real_t *uu1 = u1 + VLEN * 4 * j * ncol;
      real_t *uu2 = u2 + VLEN * 4 * j * ncol;

      // work area for field shifting
      real_t u[VLEN * 4];

      // vector variables: treat shifted varaibles
      Vsimd_t vu[4], vin[4];

      // shifted input vector
      load_vec(vin, &in[VLEN * 4 * j], 4);

      // shifted gauge field: for ch -*
      shift_vec2_yfw(u, &uu1[VLEN * 4 * i], &uu2[VLEN * 4 * i], 4);
      load_vec(vu, u, 4);

      //  out1 += conj(u[2*i  ]) * in[2*j];   // ch:--
      add_dot_vec(&out[e1r], &vu[e1r], &vin[e1r], 1);
      add_dot_vec(&out[e1r], &vu[e1i], &vin[e1i], 1);
      add_dot_vec(&out[e1i], &vu[e1r], &vin[e1i], 1);
      sub_dot_vec(&out[e1i], &vu[e1i], &vin[e1r], 1);

      // out2 -= conj(u[2*i+1]) * in[2*j];   // ch:-+
      sub_dot_vec(&out[e2r], &vu[e2r], &vin[e1r], 1);
      sub_dot_vec(&out[e2r], &vu[e2i], &vin[e1i], 1);
      sub_dot_vec(&out[e2i], &vu[e2r], &vin[e1i], 1);
      add_dot_vec(&out[e2i], &vu[e2i], &vin[e1r], 1);

      // shifted gauge field: for ch +*
      shift_vec2_yfw(u, &uu1[VLEN * (2 * ncol + 4 * i)], &uu2[VLEN * (2 * ncol + 4 * i)], 4);
      load_vec(vu, u, 4);

      // out1 -= conj(u[ncol+2*i  ]) * in[2*j+1];  // ch:+-
      sub_dot_vec(&out[e1r], &vu[e1r], &vin[e2r], 1);
      sub_dot_vec(&out[e1r], &vu[e1i], &vin[e2i], 1);
      sub_dot_vec(&out[e1i], &vu[e1r], &vin[e2i], 1);
      add_dot_vec(&out[e1i], &vu[e1i], &vin[e2r], 1);

      // out2 += conj(u[ncol+2*i+1]) * in[2*j+1];  // ch:++
      add_dot_vec(&out[e2r], &vu[e2r], &vin[e2r], 1);
      add_dot_vec(&out[e2r], &vu[e2i], &vin[e2i], 1);
      add_dot_vec(&out[e2i], &vu[e2r], &vin[e2i], 1);
      sub_dot_vec(&out[e2i], &vu[e2i], &vin[e2r], 1);
    }
  }


  //====================================================================
  static inline void accum_mult_udag(real_t *out,
                                     real_t *in,
                                     real_t *u0,
                                     const int ncol)
  {
    const int nh = ncol / 2;

    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_udag_i(tmp, in, u0, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  //====================================================================
  static inline void set_mult_u(real_t *out,
                                real_t *in,
                                real_t *u0,
                                const int ncol)
  {
    const int nh = ncol / 2;

    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      clear_vec(tmp, 4);
      accum_mult_u_i(tmp, in, u0, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  //====================================================================
  static inline void set_mult_udag(real_t *out,
                                   real_t *in,
                                   real_t *u,
                                   const int ncol)
  {
    const int nh = ncol / 2;

    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      clear_vec(tmp, 4);
      accum_mult_udag_i(tmp, in, u, i, ncol);
      save_vec(&out[VLEN * (4 * i)], tmp, 4);
    }
  }


  //====================================================================
  static inline void accum_buf(real_t *out, real_t *in, const int ncol)
  {
    Vsimd_t vin;
    Vsimd_t vout;
    for (int i = 0; i < 2 * ncol; i++) { // 2 for complex
      load_vec(&vin, &in[VLEN * i], 1);
      load_vec(&vout, &out[VLEN * i], 1);
      add_vec(&vout, &vin, 1);
      save_vec(&out[VLEN * i], &vout, 1);
    }
  }


  //====================================================================
  static inline void copy_buf(real_t *out, real_t *in, const int ncol)
  {
    for (int i = 0; i < VLEN * 2 * ncol; i++) { // 2 for complex
      out[i] = in[i];
    }
    //    Vsimd_t vin;
    //    for(int i=0; i< 2*ncol; i++){ // 2 for complex
    //      load_vec(&vin, &in[VLEN*i], 1);
    //      save_vec(&out[VLEN*i], &vin, 1);
    //    }
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_xp
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_xp1(real_t *buf, real_t *in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &in[VLEN * 4 * i], 4);
      save_vec1_x(&buf[VLENY * 4 * i], tmp, 0, 4);
    }
  }


  static inline void mult_coarse_xp2(real_t *out,
                                     real_t *u0,
                                     real_t *in0,
                                     real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nh = ncol / 2;
    real_t *in = work;
    shift_vec0_xbw(in, in0, 2 * ncol); // 2 for complex

    for (int i = 0; i < nh; i++) {
      Vsimd_t vin[4], buf[4];
      // merge the buffer
      shift_vec1_xbw(buf, &buf0[VLENY * 4 * i], 4);
      load_vec(vin, &in[VLEN * 4 * i], 4);
      add_vec(vin, buf, 4);
      save_vec(&in[VLEN * 4 * i], vin, 4);
    }
    accum_mult_u(out, in, u0, ncol);
  }


  static inline void mult_coarse_xpb(real_t *out,
                                     real_t *u,
                                     real_t *in1, real_t *in2,
                                     const int ncol, real_t *work)
  {
    real_t *in = work;
    shift_vec2_xbw(in, in1, in2, 2 * ncol); // 2 for complex
    accum_mult_u(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_xm
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_xm1(real_t *buf, real_t *u, real_t *in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      clear_vec(tmp, 4);
      accum_mult_udag_i(tmp, in, u, i, ncol);
      save_vec1_x(&buf[VLENY * 4 * i], tmp, VLENX - 1, 4);
    }
  }


  static inline void mult_coarse_xm2(real_t *out,
                                     real_t *u0, real_t *in0,
                                     real_t *buf0,
                                     const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t vtmp[4], buf[4];

      // multipliy udag
      clear_vec(vtmp, 4);
      accum_mult_udag_i(vtmp, in0, u0, i, ncol);

      // shift the result
      real_t tmp1[4 * VLEN], tmp2[4 * VLEN];
      save_vec(tmp1, vtmp, 4);
      shift_vec0_xfw(tmp2, tmp1, 4);
      load_vec(vtmp, tmp2, 4);

      // merge the buffer
      shift_vec1_xfw(buf, &buf0[VLENY * 4 * i], 4);
      add_vec(vtmp, buf, 4);

      // accumulate to the output
      Vsimd_t tmpout[4];
      load_vec(tmpout, &out[VLEN * 4 * i], 4);
      add_vec(tmpout, vtmp, 4);
      save_vec(&out[VLEN * 4 * i], tmpout, 4);
    }
  }


  static inline void mult_coarse_xmb(real_t *out,
                                     real_t *u1, real_t *u2,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    real_t *in = work;
    shift_vec2_xfw(in, in1, in2, 2 * ncol); // 2 for complex

    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_udag_xm_i(tmp, in, u1, u2, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_yp
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_yp1(real_t *buf, real_t *in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &in[VLEN * 4 * i], 4);
      save_vec1_y(&buf[VLENX * 4 * i], tmp, 0, 4);
    }
  }


  static inline void mult_coarse_yp2(real_t *out,
                                     real_t *u0,
                                     real_t *in0,
                                     real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nh = ncol / 2;
    real_t *in = work;
    shift_vec0_ybw(in, in0, 2 * ncol); // 2 for complex

    for (int i = 0; i < nh; i++) {
      Vsimd_t vin[4], buf[4];
      // merge the buffer
      shift_vec1_ybw(buf, &buf0[VLENX * 4 * i], 4);
      load_vec(vin, &in[VLEN * 4 * i], 4);
      add_vec(vin, buf, 4);
      save_vec(&in[VLEN * 4 * i], vin, 4);
    }
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_u_i(tmp, in, u0, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  static inline void mult_coarse_ypb(real_t *out,
                                     real_t *u,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    real_t *in = work;
    shift_vec2_ybw(in, in1, in2, 2 * ncol); // 2 for complex
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_u_i(tmp, in, u, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_ym
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_ym1(real_t *buf, real_t *u, real_t *in, const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      clear_vec(tmp, 4);
      accum_mult_udag_i(tmp, in, u, i, ncol);
      save_vec1_y(&buf[VLENX * 4 * i], tmp, VLENY - 1, 4);
    }
  }


  static inline void mult_coarse_ym2(real_t *out,
                                     real_t *u0, real_t *in0,
                                     real_t *buf0,
                                     const int ncol)
  {
    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t vtmp[4], buf[4];

      // multipliy udag
      clear_vec(vtmp, 4);
      accum_mult_udag_i(vtmp, in0, u0, i, ncol);

      // shift the result
      real_t tmp1[4 * VLEN], tmp2[4 * VLEN];
      save_vec(tmp1, vtmp, 4);
      shift_vec0_yfw(tmp2, tmp1, 4);
      load_vec(vtmp, tmp2, 4);

      // merge the buffer
      shift_vec1_yfw(buf, &buf0[VLENX * 4 * i], 4);
      add_vec(vtmp, buf, 4);

      // accumulate to the output
      Vsimd_t tmpout[4];
      load_vec(tmpout, &out[VLEN * 4 * i], 4);
      add_vec(tmpout, vtmp, 4);
      save_vec(&out[VLEN * 4 * i], tmpout, 4);
    }
  }


  static inline void mult_coarse_ymb(real_t *out,
                                     real_t *u1, real_t *u2,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    real_t *in = work;
    shift_vec2_yfw(in, in1, in2, 2 * ncol); // 2 for complex

    const int nh = ncol / 2;
    for (int i = 0; i < nh; i++) {
      Vsimd_t tmp[4];
      load_vec(tmp, &out[VLEN * 4 * i], 4);
      accum_mult_udag_ym_i(tmp, in, u1, u2, i, ncol);
      save_vec(&out[VLEN * 4 * i], tmp, 4);
    }
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_zp
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_zp1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


  static inline void mult_coarse_zp2(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


  static inline void mult_coarse_zpb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_zm
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_zm1(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    set_mult_udag(out, buf, u, ncol);
  }


  static inline void mult_coarse_zm2(real_t *out, real_t *buf, const int ncol)
  {
    accum_buf(out, buf, ncol);
  }


  static inline void mult_coarse_zmb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_udag(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_tp
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_tp1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


  static inline void mult_coarse_tp2(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


  static inline void mult_coarse_tpb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_tm
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_tm1(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    set_mult_udag(out, buf, u, ncol);
  }


  static inline void mult_coarse_tm2(real_t *out, real_t *buf, const int ncol)
  {
    accum_buf(out, buf, ncol);
  }


  static inline void mult_coarse_tmb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_udag(out, in, u, ncol);
  }
} // anonymous namespace

//====================================================================

template<typename AFIELD>
const std::string AFopr_Clover_coarse<AFIELD>::class_name
  = "AFopr_Clover_coarse";

//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_repr = "Dirac";  // now only the Dirac repr is available.

  int req_comm = 1;  // set 1 if communication forced any time
  //int req_comm = 0;  // set 0 if communication called in necessary

  int Ndim = CommonParameters::Ndim();

  do_comm_any = 0;
  for (int mu = 0; mu < Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("  do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_bdsize.resize(Ndim);

  int fine_nvol = CommonParameters::Nvol();
  int Nc        = CommonParameters::Nc();
  int Nd        = CommonParameters::Nd();
  int NinF      = 2 * Nc * Nd;
  workvec1.reset(NinF, fine_nvol, 1);
  workvec2.reset(NinF, fine_nvol, 1);
  workvec3.reset(NinF, fine_nvol, 1);

  // rest the timers
#ifdef AFOPR_CLOVER_COARSE_TIMER
  timer_mult.reset(new            Timer("afopr_Clover_coarse: mult           "));
  timer_pack.reset(new            Timer("afopr_Clover_coarse: pack           "));
  timer_bulk.reset(new            Timer("afopr_Clover_coarse: bulk           "));
  timer_boundary.reset(new        Timer("afopr_Clover_coarse: boundary       "));
  timer_comm.reset(new            Timer("afopr_Clover_coarse: comm           "));
  timer_comm_recv_wait.reset(new  Timer("afopr_Clover_coarse: comm_recv_wait "));
  timer_comm_send_wait.reset(new  Timer("afopr_Clover_coarse: comm_send_wait "));
  timer_comm_recv_start.reset(new Timer("afopr_Clover_coarse: comm_recv_start"));
  timer_comm_send_start.reset(new Timer("afopr_Clover_coarse: comm_send_start"));
  timer_comm_test_all.reset(new   Timer("afopr_Clover_coarse: comm_test_all  "));
  timer_clear.reset(new           Timer("afopr_Clover_coarse: clear          "));
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::setup_channels()
{
  int Ndim = CommonParameters::Ndim();

  chsend_up.resize(Ndim);
  chrecv_up.resize(Ndim);
  chsend_dn.resize(Ndim);
  chrecv_dn.resize(Ndim);

  for (int mu = 0; mu < Ndim; ++mu) {
    size_t Nvsize = m_bdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::tidyup()
{
  vout.general("%s: tidyup\n", class_name.c_str());
  ThreadManager::assert_single_thread(class_name);
  int nthreads = ThreadManager::get_num_threads_available();
  for (int i = 0; i < work_shifted.size(); i++) {
    free(work_shifted[i]);
    work_shifted[i] = nullptr;
  }

#ifdef AFOPR_CLOVER_COARSE_TIMER
  timer_mult->report();
  timer_clear->report();
  timer_pack->report();
  timer_bulk->report();
  timer_boundary->report();
  timer_comm->report();
  timer_comm_recv_wait->report();
  timer_comm_send_wait->report();
  timer_comm_recv_start->report();
  timer_comm_send_start->report();
  timer_comm_test_all->report();
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int              num_testvectors;
  std::vector<int> coarse_lattice;

  int err = 0;
  err += params.fetch_int("number_of_testvectors", num_testvectors);
  err += params.fetch_int_vector("coarse_lattice_size", coarse_lattice);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(num_testvectors, coarse_lattice);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::set_parameters(
  const int num_testvectors,
  const std::vector<int>& coarse_lattice)
{
  ThreadManager::assert_single_thread(class_name);

  int Ndim = CommonParameters::Ndim();
  assert(coarse_lattice.size() == Ndim);

  m_num_testvectors = num_testvectors;
  m_ncol            = 2 * num_testvectors; // number of chirality is multplied.
  m_Nc  = m_ncol;
  m_Nc2 = m_ncol * m_ncol;
  m_Nvc = 2 * m_Nc;        // 2 for complex
  m_Ndf = 2 * m_Nc * m_Nc; // 2 for complex
  int Nc2 = m_ncol * m_ncol;

  m_Nx   = coarse_lattice[0];
  m_Ny   = coarse_lattice[1];
  m_Nz   = coarse_lattice[2];
  m_Nt   = coarse_lattice[3];
  m_Nst  = m_Nx * m_Ny * m_Nz * m_Nt;
  m_Nstv = m_Nst / VLEN;
  m_Nxv  = m_Nx / VLENX;
  m_Nyv  = m_Ny / VLENY;

  // sanity check
  if (m_Nxv * VLENX != m_Nx) {
    vout.crucial("%s: bad coarse lattice size in x-direction: must be a multiple of %d (given: %d)\n", class_name.c_str(), VLENX, m_Nx);
    exit(EXIT_FAILURE);
  }
  if (m_Nyv * VLENY != m_Ny) {
    vout.crucial("%s: bad coarse lattice size in y-direction: must be a multiple of %d (given: %d)\n", class_name.c_str(), VLENY, m_Ny);
    exit(EXIT_FAILURE);
  }


  m_bdsize[0] = m_Nvc * m_Ny * m_Nz * m_Nt;
  m_bdsize[1] = m_Nvc * m_Nx * m_Nz * m_Nt;
  m_bdsize[2] = m_Nvc * m_Nx * m_Ny * m_Nt;
  m_bdsize[3] = m_Nvc * m_Nx * m_Ny * m_Nz;

  setup_channels();

  size_t coarse_nvol = m_Nst;
  m_coarse_lvol = coarse_nvol * CommonParameters::NPE();

  m_U.reset(m_Ndf, m_Nst, Ndim);   // hopping term
  m_Clov.reset(m_Ndf, m_Nst, 1);   // on-site term

  tmp_buffer1.resize(coarse_nvol);
  tmp_buffer2.resize(coarse_nvol);

  int nthreads  = ThreadManager::get_num_threads_available();
  int pool_size = ((sizeof(real_t) * VLEN * m_Nvc - 1) / alignment + 1) * alignment;

  for (int i = 0; i < work_shifted.size(); i++) {
    free(work_shifted[i]);
    work_shifted[i] = nullptr;
  }
  work_shifted.resize(nthreads);
  for (int i = 0; i < nthreads; i++) {
    posix_memalign((void **)&work_shifted[i], alignment, pool_size);
  }
  vout.detailed(m_vl, "shifted buffer: size=%d, alignment=%d\n", pool_size, alignment);

  set_list();
  vout.detailed(m_vl, "setting list vector, done\n");
  for (int th = 0; th < m_list_boundary.size(); th++) {
    vout.detailed(m_vl, "  thread=%d,  number of boundary sites = %d\n", th, m_list_boundary[th].size());
  }
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  coarse_lattice_size[%d] = %2d\n",
                 mu, coarse_lattice[mu]);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::set_list()
{
  int              work_xp = m_ncol;
  int              work_xm = m_ncol;
  int              work_yp = m_ncol;
  int              work_ym = m_ncol;
  int              work_zp = m_ncol;
  int              work_zm = 1;
  int              work_tp = m_ncol;
  int              work_tm = 1;
  std::vector<int> workload(m_Nstv);
  for (int site = 0; site < m_Nstv; ++site) {
    workload[site] = 0;
  }
  for (int site = 0; site < m_Nstv; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;
    int iy   = iyzt % m_Nyv;
    int izt  = site / (m_Nxv * m_Nyv);
    int iz   = izt % m_Nz;
    int it   = izt / m_Nz;
    if (do_comm[0] == 1) {
      if (ix == m_Nxv - 1) { workload[site] += work_xp; }
      if (ix == 0) { workload[site] += work_xm; }
    } // do_comm[0] == 1
    if (do_comm[1] == 1) {
      if (iy == m_Nyv - 1) { workload[site] += work_yp; }
      if (iy == 0) { workload[site] += work_ym; }
    } // do_comm[1] == 1
    if (do_comm[2] == 1) {
      if (iz == m_Nz - 1) { workload[site] += work_zp; }
      if (iz == 0) { workload[site] += work_zm; }
    } // do_comm[2] == 1
    if (do_comm[3] == 1) {
      if (it == m_Nt - 1) { workload[site] += work_tp; }
      if (it == 0) { workload[site] += work_tm; }
    } // do_comm[3] == 1
  }   // site

  int nth0 = ThreadManager::get_num_threads_available();
  int nth  = nth0;
  if (nth > 2) { nth--; }  // do not use master thread

  std::vector<std::vector<int> > tmp_list;
  std::vector<int>               work(nth);
  std::vector<int>               tmp_list_next(nth);
  tmp_list.resize(nth);
  for (int i = 0; i < nth; i++) {
    tmp_list[i].resize(m_Nstv);
  }
  for (int i = 0; i < nth; i++) {
    work[i] = 0;
  }
  for (int i = 0; i < nth; i++) {
    tmp_list_next[i] = 0;
  }
  int th_min_work = 0;
  while (1)
  {
    // find the next site, which has the maximum laod
    int max_work      = 0;
    int max_work_site = -1;
    for (int site = 0; site < m_Nstv; site++) {
      if (workload[site] > max_work) {
        max_work      = workload[site];
        max_work_site = site;
      }
    }
    if (max_work == 0) {  // no work is left
      break;
    }
    // assign the work a thread with the minnum works so far
    tmp_list[th_min_work][tmp_list_next[th_min_work]] = max_work_site;
    tmp_list_next[th_min_work]++;
    work[th_min_work]      += max_work;
    workload[max_work_site] = 0;

    // look for the next thread to work
    int min_work = work[th_min_work];
    for (int th = 0; th < nth; th++) {
      if (work[th] < min_work) {
        min_work    = work[th];
        th_min_work = th;
      }
    }
  }

  // resize and set the list vector
  m_list_boundary.resize(nth0);
  m_list_boundary[0].resize(0);
  for (int th = 0; th < nth; th++) {
    int th0 = th;
    if (nth0 > nth) { th0++; }
    int size = tmp_list_next[th];
    m_list_boundary[th0].resize(size);
    vout.general("hoge: setting list: th0=%d/%d, size=%d, load=%d\n",
                 th0, nth0, size, work[th]);
    for (int i = 0; i < size; i++) {
      vout.general("                    th0=%d/%d, i=%d site=%d\n",
                   th0, nth0, i, tmp_list[th][i]);

      m_list_boundary[th0][i] = tmp_list[th][i];
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::generate_coarse_op(
  AFopr_dd<AFIELD> *fine_afopr_,
  const std::vector<AFIELD>& atestvec)
{
  int       ith, nth, coarse_is, coarse_ns;
  const int coarse_nvol = m_U.nvol();
  set_threadtask(ith, nth, coarse_is, coarse_ns, coarse_nvol);

  real_t    *out_clov   = m_Clov.ptr(0);
  real_t    *out_gauge  = m_U.ptr(0);
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;
  assert(m_Nc2 == coarse_Nc2);

  // must be AFopr_Clover_dd
  AFopr_Clover_dd<AFIELD> *fine_afopr
    = dynamic_cast<AFopr_Clover_dd<AFIELD> *>(fine_afopr_);
  if (fine_afopr == nullptr) {
    vout.crucial("%s: in generate_coarse_op, a bad fine operator"
                 " is given (mustbbe AFopr_Clover_dd<AFIELD>).\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Clov.set(0.0);
  m_U.set(0.0);

  std::vector<int> coarse_lattice(4);
  coarse_lattice[0] = m_Nx;
  coarse_lattice[1] = m_Ny;
  coarse_lattice[2] = m_Nz;
  coarse_lattice[3] = m_Nt;
  AIndex_block_lex<real_t, QXS> index_block(coarse_lattice);

  AIndex_coarse_lex<real_t, QXS> index_coarse(m_Nx, m_Ny, m_Nz, m_Nt, num_vectors, 2);

#pragma omp barrier

  // coarse clover:
  //  I = 2*i1 + chirality1
  //  J = 2*i2 + chirality2
  //  JI = (2*num_vectors)*J + I
  //  <J|D|I>

  for (int i1 = 0; i1 < num_vectors; ++i1) {
    for (int ch1 = -1; ch1 < 2; ch1 += 2) { // ch=-1,+1
      fine_afopr->project_chiral(workvec1, atestvec[i1], ch1);

      // diag block: "clover"
      fine_afopr->mult_dd(workvec2, workvec1);
      real_t *out = out_clov;
      int    I    = 2 * i1 + (ch1 + 1) / 2;

      for (int i2 = 0; i2 < num_vectors; ++i2) {
        fine_afopr->project_chiral(workvec3, workvec2, -1);
#pragma omp barrier
        block_dotc(&tmp_buffer1[0], atestvec[i2], workvec3,
                   index_block);
#pragma omp barrier
        fine_afopr->project_chiral(workvec3, workvec2, 1);
#pragma omp barrier
        block_dotc(&tmp_buffer2[0], atestvec[i2], workvec3,
                   index_block);

        int J = 2 * i2;
        for (int s = coarse_is; s < coarse_ns; s++) {
          int idx_r = index_coarse.idx_Gr(J, I, s, 0);
          int idx_i = index_coarse.idx_Gi(J, I, s, 0);
          out[idx_r] += real(tmp_buffer1[s]);
          out[idx_i] += imag(tmp_buffer1[s]);
        }

        ++J;
        for (int s = coarse_is; s < coarse_ns; s++) {
          int idx_r = index_coarse.idx_Gr(J, I, s, 0);
          int idx_i = index_coarse.idx_Gi(J, I, s, 0);
          out[idx_r] += real(tmp_buffer2[s]);
          out[idx_i] += imag(tmp_buffer2[s]);
        }
#pragma omp barrier
      } // i2

      // hopping block: "gauge"
      for (int mu = 0; mu < 4; mu++) {
        workvec2.set(0.0);
        fine_afopr->mult_dup(workvec2, workvec1, mu);
        real_t *out = out_gauge + mu * 2 * coarse_Nc2 * coarse_nvol;
        // mu comes last

        int I = 2 * i1 + (ch1 + 1) / 2;
        for (int i2 = 0; i2 < num_vectors; ++i2) {
          fine_afopr->project_chiral(workvec3, workvec2, -1);
#pragma omp barrier
          block_dotc(&tmp_buffer1[0], atestvec[i2], workvec3,
                     index_block);
#pragma omp barrier
          fine_afopr->project_chiral(workvec3, workvec2, 1);
#pragma omp barrier
          block_dotc(&tmp_buffer2[0], atestvec[i2], workvec3,
                     index_block);

          int J = 2 * i2;
          for (int s = coarse_is; s < coarse_ns; ++s) {
            int idx_r = index_coarse.idx_Gr(J, I, s, 0);
            int idx_i = index_coarse.idx_Gi(J, I, s, 0);
            out[idx_r] += real(tmp_buffer1[s]);
            out[idx_i] += imag(tmp_buffer1[s]);
          }

          ++J;
          for (int s = coarse_is; s < coarse_ns; ++s) {
            int idx_r = index_coarse.idx_Gr(J, I, s, 0);
            int idx_i = index_coarse.idx_Gi(J, I, s, 0);
            out[idx_r] += real(tmp_buffer2[s]);
            out[idx_i] += imag(tmp_buffer2[s]);
          }
#pragma omp barrier
        } // i2
      }   // mu
    }
  }       // i1, ch1
  // rescale operator as project chiral does not have 1/2
  m_U.scal(0.25);
  m_Clov.scal(0.25);

  {
    double clv2 = m_Clov.norm2();
    double u2   = m_U.norm2();
    vout.general("%s: |m_Clov|^2 = %23.15e\n", class_name.c_str(), clv2);
    vout.general("%s: |m_U|^2    = %23.15e\n", class_name.c_str(), u2);

#ifdef DEBUG
    for (int i = 0; i < 2 * num_vectors; ++i) {
      for (int j = 0; j < 2 * num_vectors; ++j) {
        int s     = 0;
        int mu    = 3;
        int idx_r = index_coarse.idx_Gr(j, i, s, mu);
        int idx_i = index_coarse.idx_Gi(j, i, s, mu);
        vout.general("i = %d  j = %d  %f  %f\n", i, j,
                     m_U.cmp(idx_r), m_U.cmp(idx_i));
        //m_Clov.cmp(idx_r), m_Clov.cmp(idx_i));
      }
    }
#endif
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::set_config(Field *u)
{
  vout.crucial(m_vl, "%s: set_config is called\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::convert(AFIELD& v, const Field& w)
{
  // no need
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  // no need
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_up(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xp(vp, wp);
  } else if (mu == 1) {
    mult_yp(vp, wp);
  } else if (mu == 2) {
    mult_zp(vp, wp);
  } else if (mu == 3) {
    mult_tp(vp, wp);
  } else {
    vout.crucial(m_vl, "%s: mult_up for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_dn(int mu, AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  if (mu == 0) {
    mult_xm(vp, wp);
  } else if (mu == 1) {
    mult_ym(vp, wp);
  } else if (mu == 2) {
    mult_zm(vp, wp);
  } else if (mu == 3) {
    mult_tm(vp, wp);
  } else {
    vout.crucial(m_vl, "%s: mult_dn for %d direction is undefined.",
                 class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Clover_coarse<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return D(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return Ddag(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    return Ddag(v, w);
  } else if (m_mode == "DdagD") {
    return DdagD(v, w);
  } else if (m_mode == "Ddag") {
    return D(v, w);
  } else if (m_mode == "H") {
    return H(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  mult_D(v, w);
  //mult_D_alt(v, w);
  //  mult_D_alt_keep(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  mult_gm5(v, w);
  D(m_v2, v);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  mult_gm5(vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_csw(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  mult_csw(vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);
  for (int s = is; s < ns; ++s) {
    real_t *out = v + s * 2 * VLEN * m_ncol;
    const real_t *in = w + s * 2 * VLEN * m_ncol;
    for (int i = 0; i < m_ncol / 2; ++i) {
      Vsimd_t tmp[4];               // 2 for chirality, 2 for complex
      load_vec(tmp, &in[4 * i * VLEN], 4);
      scal_vec(tmp, real_t(-1), 2); // ( -re, -im, +re, +im)
      save_vec(&out[4 * i * VLEN], tmp, 4);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_csw(real_t *v2, real_t *v1)
{                               // Dirac representation is assumed.
  real_t *u = m_Clov.ptr(0);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);

#pragma omp barrier
  int nv  = VLEN * m_Nvc;
  int nv2 = VLEN * m_Ndf;

  for (int site = is; site < ns; ++site) {
    accum_mult_u(&v2[nv * site], &v1[nv * site],
                 &u[nv2 * site], m_Nc);
  }

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
#pragma omp master
  {
    TIMER_mult_start;
  }

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nstv);
  const bool time_keeper = (ith == nth - 1);
  real_t     *vp         = v.ptr(0);
  real_t     *wp         = const_cast<AFIELD *>(&w)->ptr(0);
  real_t     *up         = m_U.ptr(0);
  real_t     *cp         = m_Clov.ptr(0);

  int Nsize[4] = { m_Nxv, m_Nyv, m_Nz, m_Nt };

  if (time_keeper) {
    TIMER_pack_start;
  }
  if (do_comm_any > 0) {
#pragma omp master
    {
      TIMER_comm_recv_start_start;
      chset_recv.start();
      TIMER_comm_recv_start_stop;
    }

    real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
    real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
    real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
    real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
    real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
    real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
    real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
    real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

    BridgeQXS::mult_coarse_1(buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                             buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                             up, wp, Nsize, m_ncol, do_comm);
  }
  if (time_keeper) {
    TIMER_pack_stop;
  }

  // clear(vp);  // redundant, to be deleted

#pragma omp barrier

#pragma omp master
  {
    TIMER_comm_start;
    TIMER_comm_send_start_start;
    chset_send.start();
    TIMER_comm_send_start_stop;
  }


  if (time_keeper) { TIMER_bulk_start; }
  BridgeQXS::mult_coarse_b(vp, up, cp, wp, Nsize, m_ncol, do_comm, work_shifted[ith]);
  if (time_keeper) { TIMER_bulk_stop; } // due to load imbalacne, this timer is not accuate

#pragma omp master
  {
    TIMER_comm_recv_wait_start;
    chset_recv.wait();
    TIMER_comm_recv_wait_stop;
  }
#pragma omp barrier

#pragma omp master
  {
    TIMER_comm_stop;
  }
  if (time_keeper) { TIMER_boundary_start; }
  real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
  real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
  real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
  real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
  real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
  real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
  real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
  real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

  BridgeQXS::mult_coarse_2(vp, up, wp,
                           buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                           buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                           Nsize, m_ncol, do_comm, work_shifted[ith],
                           m_list_boundary[ith]);

#pragma omp master
  {
    TIMER_comm_send_wait_start;
    chset_send.wait();
    TIMER_comm_send_wait_stop;
  }
#pragma omp barrier
  if (time_keeper) { TIMER_boundary_stop; }

#pragma omp master
  {
    TIMER_mult_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_D_alt(AFIELD& v, const AFIELD& w)
{
#pragma omp master
  {
    TIMER_mult_start;
  }

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  clear(vp);
  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);
  mult_csw(vp, wp);

#pragma omp master
  {
    TIMER_mult_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::clear(real_t *v)
{
  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  Vsimd_t vzero[2];  // 2 for complex
  clear_vec(vzero, 2);
  for (int site = is; site < ns; ++site) {
    real_t *out = v + VLEN * m_Nvc * site;
    for (int ic = 0; ic < m_Nvc; ic += 2) {
      save_vec(&out[VLEN * ic], vzero, 2);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_xp(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[0].ptr();
  real_t *buf2 = (real_t *)chrecv_up[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  // work area for shifed vectors
  real_t *work = work_shifted[ith];

#pragma omp barrier
  if (do_comm[0] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == 0) {
        int ibf = VLENY * m_Nvc * iyzt;
        mult_coarse_xp1(&buf1[ibf], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_up[0].start();
      chsend_dn[0].start();
      chrecv_up[0].wait();
      chsend_dn[0].wait();
    }
#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  } // if(do_comm[0] == 1)

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    if ((ix < m_Nxv - 1) || (do_comm[0] == 0)) {
      int nei = (ix + 1) + m_Nxv * iyzt;
      if (ix == m_Nxv - 1) nei = 0 + m_Nxv * iyzt;
      mult_coarse_xpb(&v2[VLEN * m_Nvc * site], &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site], &v1[VLEN * m_Nvc * nei], m_Nc, work);
    } else {
      int ibf = VLENY * m_Nvc * iyzt;
      mult_coarse_xp2(&v2[VLEN * m_Nvc * site], &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site], &buf2[ibf], m_Nc, work);
    }
  }

#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  // work area for shifed vectors
  real_t *work = work_shifted[ith];

#pragma omp barrier

  if (do_comm[0] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iyzt = site / m_Nxv;
      if (ix == m_Nxv - 1) {
        int ibf = VLENY * m_Nvc * iyzt;
        mult_coarse_xm1(&buf1[ibf], &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_dn[0].start();
      chsend_up[0].start();
      chrecv_dn[0].wait();
      chsend_up[0].wait();
    }
#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  } // end of if(do_comm[0] > 0)

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iyzt = site / m_Nxv;

    if ((ix > 0) || (do_comm[0] == 0)) {
      int ix2 = (ix - 1 + m_Nxv) % m_Nxv;
      int nei = ix2 + m_Nxv * iyzt;
      mult_coarse_xmb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &u[VLEN * m_Ndf * nei],
                      &v1[VLEN * m_Nvc * site], &v1[VLEN * m_Nvc * nei],
                      m_Nc, work);
    } else {
      int ibf = VLENY * m_Nvc * iyzt;
      mult_coarse_xm2(&v2[VLEN * m_Nvc * site], &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site], &buf2[ibf], m_Nc);
    }
  }

#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_yp(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[1].ptr();
  real_t *buf2 = (real_t *)chrecv_up[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier
  if (do_comm[1] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iy   = (site / m_Nxv) % m_Nyv;
      int izt  = site / Nxyv;
      int ixzt = ix + m_Nxv * izt;
      if (iy == 0) {
        int ibf = VLENX * m_Nvc * ixzt;
        mult_coarse_yp1(&buf1[ibf], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_up[1].start();
      chsend_dn[1].start();
      chrecv_up[1].wait();
      chsend_dn[1].wait();
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }  // end of if(do_comm[1] > 0)
#pragma omp master
  {
    TIMER_bulk_start;
  }

  int    thread = ThreadManager::get_thread_id();
  real_t *work  = work_shifted[thread];
  for (int site = is; site < ns; ++site) {
    int ix   = site % m_Nxv;
    int iy   = (site / m_Nxv) % m_Nyv;
    int izt  = site / Nxyv;
    int ixzt = ix + m_Nxv * izt;

    if ((iy < m_Nyv - 1) || (do_comm[1] == 0)) {
      int iy2 = (iy + 1) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      mult_coarse_ypb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site], &v1[VLEN * m_Nvc * nei],
                      m_Nc, work);
    } else {
      int ibf = VLENX * m_Nvc * ixzt;
      mult_coarse_yp2(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site], &buf2[ibf],
                      m_Nc, work);
    }
  }


#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);


#pragma omp barrier
  if (do_comm[1] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ix  = site % m_Nxv;
      int iy  = (site / m_Nxv) % m_Nyv;
      int izt = site / Nxyv;
      if (iy == m_Nyv - 1) {
        int ibf = VLENX * m_Nvc * (ix + m_Nxv * izt);
        mult_coarse_ym1(&buf1[ibf], &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_dn[1].start();
      chsend_up[1].start();
      chrecv_dn[1].wait();
      chsend_up[1].wait();
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }

#pragma omp master
  {
    TIMER_bulk_start;
  }

  int    thread = ThreadManager::get_thread_id();
  real_t *work  = work_shifted[thread];

  for (int site = is; site < ns; ++site) {
    int ix  = site % m_Nxv;
    int iy  = (site / m_Nxv) % m_Nyv;
    int izt = site / Nxyv;

    if ((iy != 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + m_Nyv) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      mult_coarse_ymb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &u[VLEN * m_Ndf * nei],
                      &v1[VLEN * m_Nvc * site], &v1[VLEN * m_Nvc * nei],
                      m_Nc, work);
    } else {
      int ibf = VLENX * m_Nvc * (ix + m_Nxv * izt);
      mult_coarse_ym2(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site],
                      &v1[VLEN * m_Nvc * site],
                      &buf2[ibf], m_Nc);
    }
  }
#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_zp(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[2].ptr();
  real_t *buf2 = (real_t *)chrecv_up[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier
  if (do_comm[2] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ixy  = site % Nxyv;
      int iz   = (site / Nxyv) % m_Nz;
      int it   = site / (Nxyv * m_Nz);
      int ixyt = ixy + Nxyv * it;
      if (iz == 0) {
        mult_coarse_zp1(&buf1[VLEN * m_Nvc * ixyt], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_up[2].start();
      chsend_dn[2].start();
      chrecv_up[2].wait();
      chsend_dn[2].wait();
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxyv;
    int iz  = (site / Nxyv) % m_Nz;
    int it  = site / (Nxyv * m_Nz);

    if ((iz != m_Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxyv * (iz2 + m_Nz * it);
      mult_coarse_zpb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * nei], m_Nc);
    } else {
      int ixyt = ixy + Nxyv * it;
      mult_coarse_zp2(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &buf2[VLEN * m_Nvc * ixyt], m_Nc);
    }
  }

#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[2] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ixy = site % Nxyv;
      int iz  = (site / Nxyv) % m_Nz;
      int it  = site / (Nxyv * m_Nz);
      if (iz == m_Nz - 1) {
        int ixyt = ixy + Nxyv * it;
        mult_coarse_zm1(&buf1[VLEN * m_Nvc * ixyt],
                        &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_dn[2].start();
      chsend_up[2].start();
      chrecv_dn[2].wait();
      chsend_up[2].wait();
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ixy = site % Nxyv;
    int iz  = (site / Nxyv) % m_Nz;
    int it  = site / (Nxyv * m_Nz);

    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxyv * (iz2 + m_Nz * it);
      mult_coarse_zmb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * nei], &v1[VLEN * m_Nvc * nei], m_Nc);
    } else {
      int ixyt = ixy + Nxyv * it;
      mult_coarse_zm2(&v2[VLEN * m_Nvc * site], &buf2[VLEN * m_Nvc * ixyt], m_Nc);
    }
  }
#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_tp(real_t *v2, real_t *v1)
{
  int idir  = 3;
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_dn[3].ptr();
  real_t *buf2 = (real_t *)chrecv_up[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);


#pragma omp barrier

  if (do_comm[3] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyzv;
      int it   = site / Nxyzv;
      if (it == 0) {
        mult_coarse_tp1(&buf1[VLEN * m_Nvc * ixyz], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_up[3].start();
      chsend_dn[3].start();
      chrecv_up[3].wait();
      chsend_dn[3].wait();
    }

#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyzv;
    int it   = site / Nxyzv;

    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyzv * it2;
      mult_coarse_tpb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * nei], m_Nc);
    } else {
      mult_coarse_tp2(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * site], &buf2[VLEN * m_Nvc * ixyz], m_Nc);
    }
  }

#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_tm1(real_t *v2, real_t *v1)
{
  int idir  = 3;
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyzv;
      int it   = site / Nxyzv;
      if (it == m_Nt - 1) {
        mult_coarse_tm1(&buf1[VLEN * m_Nvc * ixyz],
                        &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_tmb2(real_t *v2, real_t *v1)
{
  int idir  = 3;
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyzv;
    int it   = site / Nxyzv;

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyzv * it2;
      mult_coarse_tmb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * nei], &v1[VLEN * m_Nvc * nei], m_Nc);
    } else {
      mult_coarse_tm2(&v2[VLEN * m_Nvc * site], &buf2[VLEN * m_Nvc * ixyz], m_Nc);
    }
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Clover_coarse<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir  = 3;
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * idir);

#pragma omp barrier

  if (do_comm[3] > 0) {
#pragma omp master
    {
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyzv;
      int it   = site / Nxyzv;
      if (it == m_Nt - 1) {
        mult_coarse_tm1(&buf1[VLEN * m_Nvc * ixyz],
                        &u[VLEN * m_Ndf * site], &v1[VLEN * m_Nvc * site], m_Nc);
      }
    }

#pragma omp barrier

#pragma omp master
    {
      TIMER_pack_stop;
      TIMER_comm_start;
      chrecv_dn[3].start();
      chsend_up[3].start();
      chrecv_dn[3].wait();
      chsend_up[3].wait();
    }
#pragma omp barrier
#pragma omp master
    {
      TIMER_comm_stop;
    }
  }

#pragma omp master
  {
    TIMER_bulk_start;
  }

  for (int site = is; site < ns; ++site) {
    int ixyz = site % Nxyzv;
    int it   = site / Nxyzv;

    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyzv * it2;
      mult_coarse_tmb(&v2[VLEN * m_Nvc * site],
                      &u[VLEN * m_Ndf * nei], &v1[VLEN * m_Nvc * nei], m_Nc);
    } else {
      mult_coarse_tm2(&v2[VLEN * m_Nvc * site], &buf2[VLEN * m_Nvc * ixyz], m_Nc);
    }
  }

#pragma omp barrier
#pragma omp master
  {
    TIMER_bulk_stop;
  }
}


//====================================================================
template<typename AFIELD>
double AFopr_Clover_coarse<AFIELD>::flop_count(const std::string mode)
{
  // The following counting explicitly depends on the implementation.
  // It will be recalculated when the code is modified.
  // The present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

  int    Lvol = m_Nst * CommonParameters::NPE();
  double flop_site, flop;

  //  flop_site = static_cast<double>(m_Nc * m_Nc * (4 * m_Nc));
  // each of matrix mult takes 8 N^2 flops
  //   there is a room to improve this by using a property of Hermite matrix,
  //   but is not implemented yet.
  //   [28 Mar 2021 I.Kanamori]
  flop_site = static_cast<double>(5 * 8 * m_Nc * m_Nc);

  flop = flop_site * static_cast<double>(Lvol);
  if ((mode == "DdagD") || (mode == "DDdag")) flop *= 2.0;

  return flop;
}


//============================================================END=====
