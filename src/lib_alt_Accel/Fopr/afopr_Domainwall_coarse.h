/*!
      @file    afopr_Domainwall_coarse.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2587 $
*/

#ifndef ACCEL_AFOPR_DOMAINWALL_COARSE_INCLUDED
#define ACCEL_AFOPR_DOMAINWALL_COARSE_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib_alt/Fopr/afopr_dd.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/Tools/timer.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Solver/asolver.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"

class Field;

template<typename AFIELD>
class AFopr_Domainwall_coarse : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;
  static const std::string class_name;

 protected:
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nstv;

  int m_num_testvectors;
  int m_ncol;
  int m_Nc, m_Nc2;
  int m_Nvc, m_Ndf;  //!< caution! not color, but just inner dof.

  int m_fine_fopr_nin, m_fine_fopr_nvol, m_fine_fopr_nex;

  size_t m_coarse_lvol;

  AFIELD m_U;
  AFIELD m_Clov;
  AFIELD m_Clov_inv;
  AFIELD m_U_unprec;

  AFIELD m_v1;

  std::string  m_repr;         //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;   //!< verbose level

  // AFopr_Domainwall_dd<AFIELD> *m_fine_afopr;

  AFIELD workvec1, workvec2, workvec3, workvec4;

  std::string  m_mode;  //!< mult mode

  AFIELD m_v2;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];
  int m_bc2[4];

  std::vector<int> m_bdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  std::vector<complex_t> tmp_buffer1;
  std::vector<complex_t> tmp_buffer2;
  std::vector<complex_t> tmp_buffer3;
  std::vector<real_t> tmp_buffer_real1;

public:
  //! constructor.
  AFopr_Domainwall_coarse(const Parameters& params) : AFopr<AFIELD>()
    { init();
      set_parameters(params); }

  AFopr_Domainwall_coarse() : AFopr<AFIELD>()
    { init(); }

  //! destructor.
  ~AFopr_Domainwall_coarse(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const int num_testvectors,
                      const std::vector<int> &coarse_lattice,
                      const int fine_fopr_nin,
                      const int fine_fopr_nvol,
                      const int fine_fopr_nex);

  //! setting gauge configuration.
  void set_config(Field* u);

  //! QXS2 version requires convert of spinor field.
  bool needs_convert(){return false; }

  //! setting mult mode.
  void set_mode(std::string mode){ m_mode = mode; }

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&, const string mode);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  void generate_coarse_op(AFopr_dd<AFIELD>*,
                          ASolver<AFIELD> *solver_PV,
                          const std::vector<AFIELD> &testvec);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count(){ return flop_count(m_mode); }

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

  void mult_csw(AFIELD&, const AFIELD&);

  void mult_csw_alt(AFIELD&, const AFIELD&);

  void mult_prec(AFIELD&, const AFIELD&);

  void mult_prec_alt(AFIELD&, const AFIELD&);


 private:
  //! initial setup.
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! inpose the boundary condition to link variable.
  //  void set_boundary();

  //! set_csw now assumes Dirac repr.
  //  void set_csw();

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

  void clear_host(real_t*);
  void mult_xp(real_t*, real_t*);
  void mult_xm(real_t*, real_t*);
  void mult_yp(real_t*, real_t*);
  void mult_ym(real_t*, real_t*);
  void mult_zp(real_t*, real_t*);
  void mult_zm(real_t*, real_t*);
  void mult_tp(real_t*, real_t*);
  void mult_tm(real_t*, real_t*);

  void mult_csw_host(real_t*, real_t*);

  void mult_prec_host(real_t*, real_t*);

  void mult_gm5_host(real_t*, real_t*);

  void make_preconditioned();


  //  void aypx(real_t, real_t*, real_t*);
  //  void gm5_aypx(real_t, real_t*, real_t*);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_coarse(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
                       "Domainwall_coarse", create_object_with_params);
    return init1;
  }
#endif

};

#endif  // ACCEL_AFOPR_DOMAINWALL_COARSE_INCLUDED
