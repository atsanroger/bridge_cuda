/*!
      @file    afopr_Domainwall_coarse-tmpl.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2603 $
*/

#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_dd.h"
#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/aindex_coarse_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
#include "lib/Tools/decompose_LU_Cmplx.h"
#include "omp.h"

#define AFOPR_DOMAINWALL_COARSE_TIMER

#ifdef AFOPR_DOMAINWALL_COARSE_TIMER
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
    const int nin = ncol * 2; // 2 for complex
    for (int i = 0; i < nin; i++) {
      Vsimd_t tmp;
      load_vec(&tmp, &in[VLEN * i], 1);
      save_vec1_x(&buf[VLENY * i], &tmp, 0, 1);
    }
  }


  static inline void mult_coarse_xp2(real_t *out,
                                     real_t *u0,
                                     real_t *in0,
                                     real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2;  // 2 for complex
    real_t *in = work;

    // in(kx,      ky, z, t) = in0(kx+1, ky, z, t)
    // in(VLENX-1, ky, z, t) = 0
    //  i.e.,
    //    in0 = {v[0],v[1],v[2],..,v[N-2], v[N-1]}
    //    in  = {v[1],v[2],v[3],..,v[N-1], 0     }  // shift + 0 fill
    shift_vec0_xbw(in, in0, nin);

    for (int i = 0; i < nin; i++) {
      Vsimd_t vin, buf;
      // merge the buffer
      //   buf(kx, ky, z, t) = 0 for kx < VLENX-1
      //   buf(VLENX-1, ky, z, t) = buf0(0 ,ky, z, t)
      //   i.e., buf = {0, 0, 0, buf0[0]}
      shift_vec1_xbw(&buf, &buf0[VLENY * i], 1);

      load_vec(&vin, &in[VLEN * i], 1);
      add_vec(&vin, &buf, 1);
      save_vec(&in[VLEN * i], &vin, 1);
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
  static inline void mult_coarse_xm1(real_t *buf, real_t *in, const int ncol)
  {
    // no need to mult udag
    const int nin = ncol * 2; // 2 for complex
    for (int i = 0; i < nin; i++) {
      Vsimd_t tmp;
      load_vec(&tmp, &in[VLEN * i], 1);
      save_vec1_x(&buf[VLENY * i], &tmp, VLENX-1, 1);
    }
  }


  static inline void mult_coarse_xm2(real_t *out,
                                     real_t *u0, real_t *in0,
                                     real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2;  // 2 for complex
    real_t *in = work;
    // in(kx, ky, z, t) = in0(kx-1, ky, z, t)
    // in(0,  ky, z, t) = 0
    //  i.e.,
    //    in0 = {v[0],v[1],v[2],..,v[N-2], v[N-1]}
    //    in  = {0,   v[0],v[1],..,v[N-3], v[N-2]}  // shift + 0 fill
    shift_vec0_xfw(in, in0, nin);

    for (int i = 0; i < nin; i++) {
      Vsimd_t vin, buf;
      // merge the buffer
      //   buf(kx, ky, z, t) = 0 for kx > 0
      //   buf(0, ky, z, t) = buf0(0 ,ky, z, t)
      //   i.e., buf = {buf0[0], 0,..,0}
      shift_vec1_xfw(&buf, &buf0[VLENY * i], 1);

      load_vec(&vin, &in[VLEN * i], 1);
      add_vec(&vin, &buf, 1);
      save_vec(&in[VLEN * i], &vin, 1);
    }
    accum_mult_u(out, in, u0, ncol);
  }


  static inline void mult_coarse_xmb(real_t *out,
                                     real_t *u,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    real_t *in = work;
    shift_vec2_xfw(in, in1, in2, 2 * ncol); // 2 for complex
    accum_mult_u(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_yp
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_yp1(real_t *buf, real_t *in, const int ncol)
  {
    const int nin = ncol * 2;
    for (int i = 0; i < nin; i++) {
      Vsimd_t tmp;
      load_vec(&tmp, &in[VLEN * i], 1);
      save_vec1_y(&buf[VLENX * i], &tmp, 0, 1);
    }
  }


  static inline void mult_coarse_yp2(real_t *out,
                                     real_t *u,
                                     real_t *in0,
                                     real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2; // 2 for complex
    real_t *in = work;
    shift_vec0_ybw(in, in0, nin);

    for (int i = 0; i < nin; i++) {
      Vsimd_t vin, buf;
      // merge the buffer
      shift_vec1_ybw(&buf, &buf0[VLENX * i], 1);
      load_vec(&vin, &in[VLEN * i], 1);
      add_vec(&vin, &buf, 1);
      save_vec(&in[VLEN * i], &vin, 1);
    }
    accum_mult_u(out, in, u, ncol);
  }


  static inline void mult_coarse_ypb(real_t *out,
                                     real_t *u,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2; // 2 for complex
    real_t *in = work;
    shift_vec2_ybw(in, in1, in2, nin);
    accum_mult_u(out, in, u, ncol);
  }


  ////////////////////////////////////////////////////////////////////////
  // for mult_ym
  ////////////////////////////////////////////////////////////////////////
  static inline void mult_coarse_ym1(real_t *buf, real_t *in, const int ncol)
  {
    const int nin = ncol * 2;
    for (int i = 0; i < nin; i++) {
      Vsimd_t tmp;
      load_vec(&tmp, &in[VLEN * i], 1);
      save_vec1_y(&buf[VLENX * i], &tmp, VLENY - 1, 1);
    }
  }


  static inline void mult_coarse_ym2(real_t *out,
                                     real_t *u,
                                     real_t *in0, real_t *buf0,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2;
    real_t *in = work;
    shift_vec0_yfw(in, in0, nin);

    for (int i = 0; i < nin; i++) {
      Vsimd_t vin, buf;
      // merge the buffer
      shift_vec1_yfw(&buf, &buf0[VLENX * i], 1);

      load_vec(&vin, &in[VLEN * i], 1);
      add_vec(&vin, &buf, 1);
      save_vec(&in[VLEN * i], &vin, 1);
    }
    accum_mult_u(out, in, u, ncol);
  }



  static inline void mult_coarse_ymb(real_t *out,
                                     real_t *u,
                                     real_t *in1, real_t *in2,
                                     const int ncol,
                                     real_t *work)
  {
    const int nin = ncol * 2;
    real_t *in = work;
    shift_vec2_yfw(in, in1, in2, 2 * ncol); // 2 for complex
    accum_mult_u(out, in, u, ncol);
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
  static inline void mult_coarse_zm1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


  static inline void mult_coarse_zm2(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


  static inline void mult_coarse_zmb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
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
  static inline void mult_coarse_tm1(real_t *out, real_t *in, const int ncol)
  {
    copy_buf(out, in, ncol);
  }


  static inline void mult_coarse_tm2(real_t *out, real_t *u, real_t *buf, const int ncol)
  {
    accum_mult_u(out, buf, u, ncol);
  }


  static inline void mult_coarse_tmb(real_t *out, real_t *u, real_t *in, const int ncol)
  {
    accum_mult_u(out, in, u, ncol);
  }
} // anonymous namespace

//====================================================================

template<typename AFIELD>
const std::string AFopr_Domainwall_coarse<AFIELD>::class_name
  = "AFopr_Domainwall_coarse";

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::init()
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

  // rest the timers
#ifdef AFOPR_DOMAINWALL_COARSE_TIMER
  timer_mult.reset(new            Timer("afopr_Domainwall_coarse: mult           "));
  timer_pack.reset(new            Timer("afopr_Domainwall_coarse: pack           "));
  timer_bulk.reset(new            Timer("afopr_Domainwall_coarse: bulk           "));
  timer_boundary.reset(new        Timer("afopr_Domainwall_coarse: boundary       "));
  timer_comm.reset(new            Timer("afopr_Domainwall_coarse: comm           "));
  timer_comm_recv_wait.reset(new  Timer("afopr_Domainwall_coarse: comm_recv_wait "));
  timer_comm_send_wait.reset(new  Timer("afopr_Domainwall_coarse: comm_send_wait "));
  timer_comm_recv_start.reset(new Timer("afopr_Domainwall_coarse: comm_recv_start"));
  timer_comm_send_start.reset(new Timer("afopr_Domainwall_coarse: comm_send_start"));
  timer_comm_test_all.reset(new   Timer("afopr_Domainwall_coarse: comm_test_all  "));
  timer_clear.reset(new           Timer("afopr_Domainwall_coarse: clear          "));
#endif
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::setup_channels()
{
  int Ndim = CommonParameters::Ndim();

  chsend_up.resize(Ndim);
  chrecv_up.resize(Ndim);
  chsend_dn.resize(Ndim);
  chrecv_dn.resize(Ndim);

  for (int mu = 0; mu < Ndim; ++mu) {
    size_t Nvsize = m_bdsize[mu] * sizeof(real_t);

    //    vout.general("hoge: mu=%d, m_bdsize=%d, Nvsize=%d\n", mu, m_bdsize[mu], Nvsize);
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
void AFopr_Domainwall_coarse<AFIELD>::tidyup()
{
  vout.general("%s: tidyup\n", class_name.c_str());
  ThreadManager::assert_single_thread(class_name);
  int nthreads = ThreadManager::get_num_threads_available();
  for (int i = 0; i < work_shifted.size(); i++) {
    free(work_shifted[i]);
    work_shifted[i] = nullptr;
  }

#ifdef AFOPR_DOMAINWALL_COARSE_TIMER
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
void AFopr_Domainwall_coarse<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int              num_testvectors;
  std::vector<int> coarse_lattice;
  int fine_fopr_nin, fine_fopr_nvol, fine_fopr_nex;

  int err = 0;
  err += params.fetch_int("number_of_testvectors", num_testvectors);
  err += params.fetch_int_vector("coarse_lattice_size", coarse_lattice);
  err += params.fetch_int("fine_operator_nin",  fine_fopr_nin);
  err += params.fetch_int("fine_operator_nvol", fine_fopr_nvol);
  err += params.fetch_int("fine_operator_nex",  fine_fopr_nex);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(num_testvectors, coarse_lattice,
                 fine_fopr_nin, fine_fopr_nvol, fine_fopr_nex);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_parameters(
  const int num_testvectors,
  const std::vector<int>& coarse_lattice,
  const int fine_nin, const int fine_nvol, const int fine_nex)
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

  // information of the fine grid operator, needed to generate working vestors
  m_fine_nin  = fine_nin;
  m_fine_nvol = fine_nvol;
  m_fine_nex  = fine_nex;

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

  m_U_unprec.reset(m_Ndf, m_Nst, Ndim*2);   // hopping term
  m_Clov.reset(m_Ndf, m_Nst, 1);   // on-site term

  m_U.reset(m_Ndf, m_Nst, Ndim*2);   // hopping term
  m_Clov_inv.reset(m_Ndf, m_Nst, 1);   // on-site term

  tmp_buffer1.resize(coarse_nvol);
  tmp_buffer2.resize(coarse_nvol);
  tmp_buffer_real1.resize(coarse_nvol);

  workvec1.reset(m_fine_nin, m_fine_nvol, m_fine_nex);
  workvec2.reset(m_fine_nin, m_fine_nvol, m_fine_nex);
  workvec3.reset(m_fine_nin, m_fine_nvol, m_fine_nex);

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

  //  set_list();
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


/*
//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_list()
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
*/

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::project_chiral(
                                          AFIELD& v, const AFIELD& w,
                                          const int ch,
                                          AFopr<AFIELD> *fine_afopr)
{
  fine_afopr->mult_gm5(v, w);
  if (ch == 1) {
    ::aypx(real_t(1.0), v, w);
  } else {
    ::aypx(real_t(-1.0), v, w);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::generate_coarse_op(
                                AFopr_dd<AFIELD> *fine_afopr_,
                                ASolver<AFIELD> *solver_PV,
                                const std::vector<AFIELD>& atestvec)
{
  int ith, nth, coarse_is, coarse_ns;
  const int coarse_nvol = m_U.nvol();
  set_threadtask(ith, nth, coarse_is, coarse_ns, coarse_nvol);

  real_t    *out_clov   = m_Clov.ptr(0);
  real_t    *out_gauge  = m_U_unprec.ptr(0);
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;
  assert(m_Nc2 == coarse_Nc2);

  // must be AFopr_Domainwall_dd
  AFopr_Domainwall_5din_dd<AFIELD> *fine_afopr
    = dynamic_cast<AFopr_Domainwall_5din_dd<AFIELD> *>(fine_afopr_);
  if (fine_afopr == nullptr) {
    vout.crucial("%s: in generate_coarse_op, a bad fine operator"
                 " is given (mustbbe AFopr_Domainwall_5din_dd<AFIELD>).\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Clov.set(0.0);
  m_U_unprec.set(0.0);


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

#pragma omp for
  for(int i=0; i<tmp_buffer_real1.size(); ++i){
    tmp_buffer_real1[i] = 0;
  }
  for (int i1 = 0; i1 < num_vectors; ++i1) {
    for (int ch1 = -1; ch1 < 2; ch1 += 2) { // ch=-1,+1

      for(int ieo=0; ieo<2; ++ieo) {
        // diag block: "clover"
        project_chiral(workvec2, atestvec[i1], ch1, fine_afopr);

        // calculate work1 := M_{PV})_{SAP} |i1, eo>
        int ieo_mask = 1-ieo;
        int Nconv;
        real_t diff;
        block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo_mask, index_block);
        //        vout.general("hoge: calling Solver_PV\n");
        //solver_PV->solve(workvec1, workvec2, Nconv, diff);
        solver_PV->solve(workvec1, workvec2, Nconv, diff, ieo);
        /*
        {
          double w1=workvec1.norm2();
          double w2=workvec2.norm2();
          vout.general("hoge: i1=%d, ch1=%d, ieo=%d, w1=%23.15e, w2=%23.15e\n", i1, ch1, ieo, w1, w2);
        }
        */
        fine_afopr->mult(workvec2, workvec1);
        block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo_mask, index_block);

        real_t *out = out_clov;
        int    I    = 2 * i1 + (ch1 + 1) / 2;

        for (int i2 = 0; i2 < num_vectors; ++i2) {
          project_chiral(workvec3, workvec2, -1, fine_afopr);
#pragma omp barrier
          block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                        ieo, index_block);
#pragma omp barrier
          project_chiral(workvec3, workvec2, 1, fine_afopr);
#pragma omp barrier
          block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                        ieo, index_block);

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

        // hopping block: "gauge "
        for (int mu = 0; mu < 8; mu++) {
          workvec2.set(0.0);
          if( mu<4 ) {
            fine_afopr->mult_dup(workvec2, workvec1, mu);
          } else {
            fine_afopr->mult_ddn(workvec2, workvec1, mu-4);
          }
          //          block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo, index_block);

          real_t *out = out_gauge + mu * 2 * coarse_Nc2 * coarse_nvol;
          // mu comes last

          int I = 2 * i1 + (ch1 + 1) / 2;
          for (int i2 = 0; i2 < num_vectors; ++i2) {
            project_chiral(workvec3, workvec2, -1, fine_afopr);
#pragma omp barrier
            block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                          1-ieo, index_block);
#pragma omp barrier
            project_chiral(workvec3, workvec2, 1, fine_afopr);
#pragma omp barrier
            block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                          1-ieo, index_block);

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
      } // ieo
    }
  }       // i1, ch1
  // rescale operator as project chiral does not have 1/2
  m_U_unprec.scal(0.25);
  m_Clov.scal(0.25);

  {
    double clv2 = m_Clov.norm2();
    double u2   = m_U_unprec.norm2();
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
                     m_U_unprec.cmp(idx_r), m_U_unprec.cmp(idx_i));
        //m_Clov.cmp(idx_r), m_Clov.cmp(idx_i));
      }
    }
#endif
  }


  make_preconditioned();


  {
    double clv2 = m_Clov_inv.norm2();
    double u2   = m_U.norm2();
    vout.general("%s: after preconditioned operator\n", class_name.c_str());
    vout.general("%s: |m_Clov_inv|^2 = %23.15e\n", class_name.c_str(), clv2);
    vout.general("%s: |m_U|^2        = %23.15e\n", class_name.c_str(), u2);
  }
}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::generate_coarse_op(
                                 AFopr_dd<AFIELD> *fine_afopr,
                                 const std::vector<AFIELD>& atestvec)
{
  int ith, nth, coarse_is, coarse_ns;
  const int coarse_nvol = m_U.nvol();
  set_threadtask(ith, nth, coarse_is, coarse_ns, coarse_nvol);

  real_t    *out_clov   = m_Clov.ptr(0);
  real_t    *out_gauge  = m_U_unprec.ptr(0);
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;
  assert(m_Nc2 == coarse_Nc2);

  m_Clov.set(0.0);
  m_U_unprec.set(0.0);

  std::vector<int> coarse_lattice(4);
  coarse_lattice[0] = m_Nx;
  coarse_lattice[1] = m_Ny;
  coarse_lattice[2] = m_Nz;
  coarse_lattice[3] = m_Nt;
  AIndex_block_lex<real_t, QXS> index_block(coarse_lattice);

  AIndex_coarse_lex<real_t, QXS> index_coarse(m_Nx, m_Ny, m_Nz, m_Nt,
                                              num_vectors, 2);

#pragma omp barrier

  // coarse clover:
  //  I = 2*i1 + chirality1
  //  J = 2*i2 + chirality2
  //  JI = (2*num_vectors)*J + I
  //  <J|D|I>

#pragma omp for
  for(int i=0; i<tmp_buffer_real1.size(); ++i){
    tmp_buffer_real1[i] = 0;
  }
  for (int i1 = 0; i1 < num_vectors; ++i1) {
    for (int ch1 = -1; ch1 < 2; ch1 += 2) { // ch=-1,+1

      for(int ieo=0; ieo<2; ++ieo) {
        // diag block: "clover"
        //project_chiral(workvec2, atestvec[i1], ch1, fine_afopr);
        project_chiral(workvec1, atestvec[i1], ch1, fine_afopr);

        // calculate work1 := M_{PV})_{SAP} |i1, eo>
        int ieo_mask = 1-ieo;
        int Nconv;
        real_t diff;
        block_scal_eo(workvec1, &tmp_buffer_real1[0], ieo_mask, index_block);

        fine_afopr->mult(workvec2, workvec1);
        block_scal_eo(workvec2, &tmp_buffer_real1[0], ieo_mask, index_block);

        real_t *out = out_clov;
        int    I    = 2 * i1 + (ch1 + 1) / 2;

        for (int i2 = 0; i2 < num_vectors; ++i2) {
          project_chiral(workvec3, workvec2, -1, fine_afopr);
#pragma omp barrier
          block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                        ieo, index_block);
#pragma omp barrier
          project_chiral(workvec3, workvec2, 1, fine_afopr);
#pragma omp barrier
          block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                        ieo, index_block);

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

        // hopping block: "gauge "
        for (int mu = 0; mu < 8; mu++) {
          workvec2.set(0.0);
          if( mu<4 ) {
            fine_afopr->mult_dup(workvec2, workvec1, mu);
          } else {
            fine_afopr->mult_ddn(workvec2, workvec1, mu-4);
          }

          real_t *out = out_gauge + mu * 2 * coarse_Nc2 * coarse_nvol;
          // mu comes last

          int I = 2 * i1 + (ch1 + 1) / 2;
          for (int i2 = 0; i2 < num_vectors; ++i2) {
            project_chiral(workvec3, workvec2, -1, fine_afopr);
#pragma omp barrier
            block_dotc_eo(&tmp_buffer1[0], atestvec[i2], workvec3,
                          1-ieo, index_block);
#pragma omp barrier
            project_chiral(workvec3, workvec2, 1, fine_afopr);
#pragma omp barrier
            block_dotc_eo(&tmp_buffer2[0], atestvec[i2], workvec3,
                          1-ieo, index_block);

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
      } // ieo
    }
  }       // i1, ch1
  // rescale operator as project chiral does not have 1/2
  m_U_unprec.scal(0.25);
  m_Clov.scal(0.25);

  {
    double clv2 = m_Clov.norm2();
    double u2   = m_U_unprec.norm2();
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
                     m_U_unprec.cmp(idx_r), m_U_unprec.cmp(idx_i));
        //m_Clov.cmp(idx_r), m_Clov.cmp(idx_i));
      }
    }
#endif
  }


  make_preconditioned();


  {
    double clv2 = m_Clov_inv.norm2();
    double u2   = m_U.norm2();
    vout.general("%s: after preconditioned operator\n", class_name.c_str());
    vout.general("%s: |m_Clov_inv|^2 = %23.15e\n", class_name.c_str(), clv2);
    vout.general("%s: |m_U|^2        = %23.15e\n", class_name.c_str(), u2);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::make_preconditioned()
{
  int ith, nth, is, ns;
  const int coarse_nvol = m_U.nvol();
  const int num_vectors = m_num_testvectors;
  const int coarse_Nc2  = m_ncol * m_ncol;

  set_threadtask(ith, nth, is, ns, coarse_nvol);

  std::vector<int> coarse_lattice(4);
  coarse_lattice[0] = m_Nx;
  coarse_lattice[1] = m_Ny;
  coarse_lattice[2] = m_Nz;
  coarse_lattice[3] = m_Nt;

  AIndex_block_lex<real_t, QXS> index_block(coarse_lattice);
  AIndex_coarse_lex<real_t, QXS> index_coarse(
                         m_Nx, m_Ny, m_Nz, m_Nt, num_vectors, 2);

  std::vector<double> mat(2*coarse_Nc2);

  const real_t  *clov  = m_Clov.ptr(0);
  const real_t  *gauge = m_U_unprec.ptr(0);
  real_t    *clov_inv  = m_Clov_inv.ptr(0);
  real_t    *gauge_out = m_U.ptr(0);

  Decompose_LU_Cmplx lu(m_ncol);
  for(int site = is; site<ns; ++site){
    for(int i=0; i<2*coarse_Nc2; ++i) {
      mat[i] = clov[index_coarse.idx_G(i, site, 0)];
    }
    lu.set_matrix(&mat[0]);

    // preconditioning matrix
    lu.get_inverse(&mat[0]);
    for(int i=0; i<2*coarse_Nc2; ++i) {
      clov_inv[index_coarse.idx_G(i, site, 0)] = mat[i];
    }

    // preconditioned matrix
    for(int mu = 0; mu < 8; ++mu) {
      for(int i=0; i<2*coarse_Nc2; ++i) {
        mat[i] = gauge[index_coarse.idx_G(i, site, mu)];
      }
      lu.mult_inverse(&mat[0]);
      for(int i=0; i<2*coarse_Nc2; ++i) {
        gauge_out[index_coarse.idx_G(i, site, mu)] = mat[i];
      }

    }
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::set_config(Field *u)
{
  vout.crucial(m_vl, "%s: set_config is called\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::convert(AFIELD& v, const Field& w)
{
  // no need
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::reverse(Field& v, const AFIELD& w)
{
  // no need
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_up(int mu, AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_coarse<AFIELD>::mult_dn(int mu, AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_coarse<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
std::string AFopr_Domainwall_coarse<AFIELD>::get_mode() const
{
  return m_mode;
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult(AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_coarse<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                           const string mode)
{
  if (mode == "prec") {
    mult_prec(v, w);
  } else if (mode == "clov") {
    mult_clov(v, w);
  } else {
    vout.crucial(m_vl, "%s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
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
void AFopr_Domainwall_coarse<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  //  mult_D(v, w);
  mult_D_alt(v, w);
  //  mult_D_alt_keep(v, w);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: DdagD is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: Ddag is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_gm5(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_gm5 is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  mult_gm5(vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_csw(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  mult_csw(vp, wp);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_prec(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier
  clear(vp);
  mult_csw_inv(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_clov(AFIELD& v, const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier
  clear(vp);
  mult_csw(vp, wp);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_gm5(real_t *v, real_t *w)
{
  vout.crucial(m_vl, "%s: mult_gm5 is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);

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
void AFopr_Domainwall_coarse<AFIELD>::mult_csw(real_t *v2, real_t *v1)
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
void AFopr_Domainwall_coarse<AFIELD>::mult_csw_inv(real_t *v2, real_t *v1)
{                               // Dirac representation is assumed.
  real_t *u = m_Clov_inv.ptr(0);

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
void AFopr_Domainwall_coarse<AFIELD>::mult_D(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: mult_D is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);
  /*
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
  */
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_D_alt(AFIELD& v, const AFIELD& w)
{

#pragma omp barrier

#pragma omp master
  {
    TIMER_mult_start;
  }
  copy(v,w);
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  mult_xp(vp, wp);
  mult_xm(vp, wp);
  mult_yp(vp, wp);
  mult_ym(vp, wp);
  mult_zp(vp, wp);
  mult_zm(vp, wp);
  mult_tp(vp, wp);
  mult_tm(vp, wp);

#pragma omp master
  {
    TIMER_mult_stop;
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
  vout.crucial(m_vl, "%s: H is not implemented\n", class_name.c_str());
  exit(EXIT_FAILURE);
  D(m_v2, w);
  mult_gm5(v, m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::clear(real_t *v)
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
void AFopr_Domainwall_coarse<AFIELD>::mult_xp(real_t *v2, real_t *v1)
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
        int iferm  = VLEN  * m_Nvc * site;
        int ibf    = VLENY * m_Nvc * iyzt;
        mult_coarse_xp1(&buf1[ibf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;

    if ((ix < m_Nxv - 1) || (do_comm[idir] == 0)) {
      int ix2 = (ix + 1) % m_Nxv;
      int nei = ix2 + m_Nxv * iyzt;
      int iferm_nei = VLEN * m_Nvc * nei;
      mult_coarse_xpb(&v2[iferm], &u[igauge],
                      &v1[iferm], &v1[iferm_nei], m_Nc, work);
    } else {
      int ibf = VLENY * m_Nvc * iyzt;
      mult_coarse_xp2(&v2[iferm], &u[igauge],
                      &v1[iferm], &buf2[ibf], m_Nc, work);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_xm(real_t *v2, real_t *v1)
{
  int idir = 0;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[0].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[0].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir+4));

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
        int iferm  = VLEN * m_Nvc * site;
        mult_coarse_xm1(&buf1[ibf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;

    if ((ix != 0) || (do_comm[idir] == 0)) {
      int ix2 = (ix - 1 + m_Nxv) % m_Nxv;
      int nei = ix2 + m_Nxv * iyzt;
      int iferm_nei  = VLEN * m_Nvc * nei;
      mult_coarse_xmb(&v2[iferm], &u[igauge],
                      &v1[iferm], &v1[iferm_nei], m_Nc, work);
    } else {
      int ibf = VLENY * m_Nvc * iyzt;
      mult_coarse_xm2(&v2[iferm], &u[igauge],
                      &v1[iferm], &buf2[ibf], m_Nc, work);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_yp(real_t *v2, real_t *v1)
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
      Communicator::sync();
      TIMER_pack_start;
    }

    for (int site = is; site < ns; ++site) {
      int ix   = site % m_Nxv;
      int iy   = (site / m_Nxv) % m_Nyv;
      int izt  = site / Nxyv;
      int ixzt = ix + m_Nxv * izt;

      if (iy == 0) {
        int iferm  = VLEN  * m_Nvc * site;
        int ibf    = VLENX * m_Nvc * ixzt;
        mult_coarse_yp1(&buf1[ibf], &v1[iferm], m_Nc);
      }
    }
#pragma omp master
    {
      Communicator::sync();
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
#pragma omp barrier
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;
    if ((iy < m_Nyv - 1) || (do_comm[1] == 0)) {
      int iy2 = (iy + 1) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      int iferm_nei  = VLEN * m_Nvc * nei;
      mult_coarse_ypb(&v2[iferm], &u[igauge],
                      &v1[iferm], &v1[iferm_nei], m_Nc, work);
    } else {
      int ibf = VLENX * m_Nvc * ixzt;
      mult_coarse_yp2(&v2[iferm], &u[igauge],
                      &v1[iferm], &buf2[ibf], m_Nc, work);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_ym(real_t *v2, real_t *v1)
{
  int idir = 1;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[1].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[1].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir+4));


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
        int ibf   = VLENX * m_Nvc * (ix + m_Nxv * izt);
        int iferm = VLEN  * m_Nvc * site;
        mult_coarse_ym1(&buf1[ibf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;

    if ((iy != 0) || (do_comm[idir] == 0)) {
      int iy2 = (iy - 1 + m_Nyv) % m_Nyv;
      int nei = ix + m_Nxv * (iy2 + m_Nyv * izt);
      int iferm_nei = VLEN * m_Nvc * nei;
      mult_coarse_ymb(&v2[iferm], &u[igauge],
                      &v1[iferm], &v1[iferm_nei], m_Nc, work);
    } else {
      int ibf = VLENX * m_Nvc * (ix + m_Nxv * izt);
      mult_coarse_ym2(&v2[iferm], &u[igauge],
                      &v1[iferm], &buf2[ibf], m_Nc, work);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_zp(real_t *v2, real_t *v1)
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
        int ibuf  = VLEN * m_Nvc * ixyt;
        int iferm = VLEN * m_Nvc * site;
        mult_coarse_zp1(&buf1[ibuf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;
    if ((iz != m_Nz - 1) || (do_comm[2] == 0)) {
      int iz2 = (iz + 1) % m_Nz;
      int nei = ixy + Nxyv * (iz2 + m_Nz * it);
      int iferm_nei = VLEN * m_Nvc * nei;
      mult_coarse_zpb(&v2[iferm], &u[igauge], &v1[iferm_nei], m_Nc);
    } else {
      int ixyt = ixy + Nxyv * it;
      int ibuf = VLEN * m_Nvc * ixyt;
      mult_coarse_zp2(&v2[iferm], &u[igauge], &buf2[ibuf], m_Nc);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_zm(real_t *v2, real_t *v1)
{
  int idir = 2;
  int Nxyv = m_Nxv * m_Nyv;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[2].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[2].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir+4));

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
        int ibuf = VLEN * m_Nvc * ixyt;
        int iferm = VLEN * m_Nvc * site;
        mult_coarse_zm1(&buf1[ibuf], &v1[iferm], m_Nc);
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
    int iferm =  VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;
    if ((iz > 0) || (do_comm[2] == 0)) {
      int iz2 = (iz - 1 + m_Nz) % m_Nz;
      int nei = ixy + Nxyv * (iz2 + m_Nz * it);
      int iferm_nei = VLEN * m_Nvc * nei;
      mult_coarse_zmb(&v2[iferm], &u[igauge], &v1[iferm_nei], m_Nc);
    } else {
      int ixyt = ixy + Nxyv * it;
      int ibuf = VLEN * m_Nvc * ixyt;
      mult_coarse_zm2(&v2[iferm], &u[igauge], &buf2[ibuf], m_Nc);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_tp(real_t *v2, real_t *v1)
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
        int ibuf =  VLEN * m_Nvc * ixyz;
        int iferm = VLEN * m_Nvc * site;
        mult_coarse_tp1(&buf1[ibuf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;
    if ((it < m_Nt - 1) || (do_comm[3] == 0)) {
      int it2 = (it + 1) % m_Nt;
      int nei = ixyz + Nxyzv * it2;
      int iferm_nei = VLEN * m_Nvc * nei;
      mult_coarse_tpb(&v2[iferm], &u[igauge], &v1[iferm_nei], m_Nc);
    } else {
      int ibuf = VLEN * m_Nvc * ixyz;
      mult_coarse_tp2(&v2[iferm], &u[igauge], &buf2[ibuf], m_Nc);
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
void AFopr_Domainwall_coarse<AFIELD>::mult_tm1(real_t *v2, real_t *v1)
{
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();

  if (do_comm[3] > 0) {
    for (int site = is; site < ns; ++site) {
      int ixyz = site % Nxyzv;
      int it   = site / Nxyzv;
      if (it == m_Nt - 1) {
        int ibuf  = VLEN * m_Nvc * ixyz;
        int iferm = VLEN * m_Nvc * site;
        mult_coarse_tm1(&buf1[VLEN * m_Nvc * ixyz], &v1[iferm], m_Nc);
      }
    }
  }
}




//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_coarse<AFIELD>::mult_tm(real_t *v2, real_t *v1)
{
  int idir  = 3;
  int Nxyzv = m_Nxv * m_Nyv * m_Nz;

  int ith, nth, is, ns;
  set_threadtask_mult(ith, nth, is, ns, m_Nstv);

  real_t *buf1 = (real_t *)chsend_up[3].ptr();
  real_t *buf2 = (real_t *)chrecv_dn[3].ptr();

  real_t *u = m_U.ptr(m_Ndf * m_Nst * (idir+4));

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
        int ibuf = VLEN * m_Nvc * ixyz;
        int iferm = VLEN * m_Nvc * site;
        mult_coarse_tm1(&buf1[ibuf], &v1[iferm], m_Nc);
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
    int iferm  = VLEN * m_Nvc * site;
    int igauge = VLEN * m_Ndf * site;
    if ((it > 0) || (do_comm[3] == 0)) {
      int it2 = (it - 1 + m_Nt) % m_Nt;
      int nei = ixyz + Nxyzv * it2;
      int iferm_nei  = VLEN * m_Nvc * nei;
      mult_coarse_tmb(&v2[iferm], &u[igauge], &v1[iferm_nei], m_Nc);
    } else {
      int ibuf = VLEN * m_Nvc * ixyz;
      mult_coarse_tm2(&v2[iferm], &u[igauge], &buf2[ibuf], m_Nc);
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
double AFopr_Domainwall_coarse<AFIELD>::flop_count(const std::string mode)
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
