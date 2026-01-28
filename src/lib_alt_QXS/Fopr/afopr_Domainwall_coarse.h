/*!
      @file    afopr_Domainwall_coarse.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2603 $
*/

#ifndef QXS_AFOPR_DOMAINWALL_COARSE_INCLUDED
#define QXS_AFOPR_DOMAINWALL_COARSE_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt/Solver/asolver.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

class Field;

template<typename AFIELD>
class AFopr_Domainwall_coarse : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t      real_t;
  typedef typename AFIELD::complex_t   complex_t;
  static const std::string class_name;

 protected:
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nstv;

  int m_num_testvectors;
  int m_ncol;
  int m_Nc, m_Nc2;
  int m_Nvc, m_Ndf;  //!< caution! not color, but just inner dof.

  int m_fine_nin;
  int m_fine_nvol;
  int m_fine_nex;

  size_t m_coarse_lvol;

  AFIELD m_U_unprec;
  AFIELD m_Clov;
  AFIELD m_Clov_inv;
  AFIELD m_U;

  AFIELD m_v1;

  std::string m_repr;         //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;  //!< verbose level

  AFIELD workvec1, workvec2, workvec3;

  std::string m_mode;   //!< mult mode

  AFIELD m_v2;


  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_bdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  std::vector<complex_t> tmp_buffer1;
  std::vector<complex_t> tmp_buffer2;
  std::vector<real_t>    tmp_buffer_real1; // real
  std::vector<real_t *> work_shifted;

 public:
  //! constructor.
  AFopr_Domainwall_coarse(const Parameters& params) : AFopr<AFIELD>()
  {
    init();
    set_parameters(params);
  }

  AFopr_Domainwall_coarse() : AFopr<AFIELD>()
  {
    init();
  }

  //! destructor.
  ~AFopr_Domainwall_coarse() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const int num_testvectors,
                      const std::vector<int>& coarse_lattice,
                      const int fine_nin, const int fine_nvol, const int fine_nex);

  //! setting gauge configuration.
  void set_config(Field *u);

  //! QXS version requires convert of spinor field.
  bool needs_convert() { return true; }

  //! convert of spinor field.
  void convert(AFIELD& v, const Field& w);

  //! reverse of spinor field.
  void reverse(Field& v, const AFIELD& w);

  //! returns the pointer to gauge configuration.
  //  inline Field* get_conf(void){return &m_U.ptr();};

  //! setting mult mode.
  //  void set_mode(std::string mode){ m_mode = mode; }
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult(AFIELD&, const AFIELD&, const string mode);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_csw(AFIELD&, const AFIELD&);
  void mult_prec(AFIELD&, const AFIELD&);
  void mult_clov(AFIELD&, const AFIELD&);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  void generate_coarse_op(AFopr_dd<AFIELD> *fine_afopr,
                          ASolver<AFIELD> *solver_PV,
                          const std::vector<AFIELD>& testvec);

  void generate_coarse_op(AFopr_dd<AFIELD> *fine_afopr,
                          const std::vector<AFIELD>& testvec);


  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count() { return flop_count(m_mode); }

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

  // for debug
  void dump();

 private:
  //! initial setup.
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! inpose the boundary condition to link variable.
  void set_boundary();

  //! prepare preconditioned operator
  void make_preconditioned();

  void project_chiral(AFIELD& v, const AFIELD& w, const int ch, AFopr<AFIELD> *afopr_fine);

  void mult_csw(real_t *v2, real_t *v1);
  void mult_csw_inv(real_t *v2, real_t *v1);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  void mult_gm4(AFIELD&, const AFIELD&);

  //! standard D mult.
  void mult_D(AFIELD&, const AFIELD&);

  //! D mult using QWS library.
  //void mult_D_qws(AFIELD&, const AFIELD&);

  //! D mult using mult_xp, etc.
  void mult_D_alt(AFIELD&, const AFIELD&);
  void mult_D_alt_keep(AFIELD&, const AFIELD&);
  void mult_D_alt_keep2(AFIELD&, const AFIELD&);

  void mult_xp(real_t *, real_t *);
  void mult_xm(real_t *, real_t *);
  void mult_yp(real_t *, real_t *);
  void mult_ym(real_t *, real_t *);
  void mult_zp(real_t *, real_t *);
  void mult_zm(real_t *, real_t *);
  void mult_tp(real_t *, real_t *);
  void mult_tm(real_t *, real_t *);

  void mult_xp1(real_t *, real_t *);
  void mult_xm1(real_t *, real_t *);
  void mult_yp1(real_t *, real_t *);
  void mult_ym1(real_t *, real_t *);
  void mult_zp1(real_t *, real_t *);
  void mult_zm1(real_t *, real_t *);
  void mult_tp1(real_t *, real_t *);
  void mult_tm1(real_t *, real_t *);

  void mult_xpb2(real_t *, real_t *);
  void mult_xmb2(real_t *, real_t *);
  void mult_ypb2(real_t *, real_t *);
  void mult_ymb2(real_t *, real_t *);
  void mult_zpb2(real_t *, real_t *);
  void mult_zmb2(real_t *, real_t *);
  void mult_tpb2(real_t *, real_t *);
  void mult_tmb2(real_t *, real_t *);

  void mult_gm5(real_t *, real_t *);

  void clear(real_t *);

  void aypx(complex_t, real_t *, real_t *);
  void gm5_aypx(complex_t, real_t *, real_t *);

  void set_list();

  std::vector<std::vector<int> > m_list_boundary;
  unique_ptr<Timer> timer_pack;
  unique_ptr<Timer> timer_bulk;
  unique_ptr<Timer> timer_boundary;
  unique_ptr<Timer> timer_comm;
  unique_ptr<Timer> timer_comm_recv_wait;
  unique_ptr<Timer> timer_comm_send_wait;
  unique_ptr<Timer> timer_comm_send_start;
  unique_ptr<Timer> timer_comm_recv_start;
  unique_ptr<Timer> timer_comm_test_all;
  unique_ptr<Timer> timer_mult;
  unique_ptr<Timer> timer_clear;

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_coarse(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Domainwall_coarse",
                                                         create_object_with_params);
    return init1;
  }
#endif
};

#endif  // AFOPR_DOMAINWALL_INCLUDED
