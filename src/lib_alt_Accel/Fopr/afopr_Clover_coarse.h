/*!
      @file    afopr_Clover_coarse.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2229 $
*/

#ifndef ACCEL_AFOPR_CLOVER_COARSE_INCLUDED
#define ACCEL_AFOPR_CLOVER_COARSE_INCLUDED

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
#include "lib/IO/bridgeIO.h"
using Bridge::vout;


#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
//#include "lib_alt_Accel/Fopr/afopr_Clover_dd.h"

class Field;

template<typename AFIELD>
class AFopr_Clover_coarse : public AFopr<AFIELD>
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

  size_t m_coarse_lvol;

  AFIELD m_U;
  AFIELD m_Clov;
  AFIELD m_v1;

  std::string  m_repr;         //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;   //!< verbose level

  // AFopr_Clover_dd<AFIELD> *m_fine_afopr;

  AFIELD workvec1, workvec2, workvec3, workvec4;

  std::string  m_mode;  //!< mult mode

  AFIELD m_v2;

  AFIELD m_T;  //!< clover term

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

public:
  //! constructor.
  AFopr_Clover_coarse(const Parameters& params) : AFopr<AFIELD>()
    { init();
      set_parameters(params); }

  AFopr_Clover_coarse() : AFopr<AFIELD>()
    { init();
    }

  //! destructor.
  ~AFopr_Clover_coarse(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const int num_testvectors,
                      const std::vector<int> &coarse_lattice);

  //! setting gauge configuration.
  void set_config(Field* u);

  //! QXS2 version requires convert of spinor field.
  bool needs_convert(){return true; }

  //! convert of spinor field.
  void convert(AFIELD& v, const Field& w);

  //! reverse of spinor field.
  void reverse(Field& v, const AFIELD& w);
    
  //! returns the pointer to gauge configuration.
  //  inline Field* get_conf(void){return &m_U.ptr();};

  //! setting mult mode.
  void set_mode(std::string mode){ m_mode = mode; }

  //! setting fine grid operator.
  //  void set_fine_op(AFopr_Clover_dd<AFIELD> *fine_afopr){
  //    m_fine_afopr=fine_afopr;
  //  }

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  void generate_coarse_op(AFopr_dd<AFIELD>*,
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

  void mult_gm5_host(real_t*, real_t*);

  //  void aypx(real_t, real_t*, real_t*);
  //  void gm5_aypx(real_t, real_t*, real_t*);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover_coarse(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Clover_coarse",
                                              create_object_with_params);
    return init1;
  }
#endif

};

#endif  // AFOPR_CLOVER_INCLUDED
