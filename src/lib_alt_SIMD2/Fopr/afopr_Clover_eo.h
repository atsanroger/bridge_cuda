/*!
      @file    afopr_Clover_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef AFOPR_CLOVER_EO_H
#define AFOPR_CLOVER_EO_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Tools/timer.h"

#if defined USE_MPI
#include "lib/Communicator/MPI/communicator_mpi.h"
#include "lib/Communicator/MPI/channel.h"
#else
#include "lib/Communicator/Single/channel.h"
#endif

#include "lib_alt_SIMD2/Fopr/afopr_CloverTerm.h"

class Field;

/*!
  This implementation is for an architecture with
  vector arithmetic units and registers.
                               [2 Feb 2016, H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Clover_eo : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd;
  int m_Ndf, m_Ncol, m_Nvc;
  int m_Nst, m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Ndim;
  int m_Nx2, m_Nst2;
  int m_Nx2v, m_Nst2v;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_cSW;               //!< clover coefficient.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::vector<int> m_Leo;  //!< Leo = 0 (even site) or 1 (odd site).

  Field *m_conf;        //!< original gauge config.
  AFIELD m_Ueo; //!< copied gauge config. with boundary conditions.

  AFIELD m_fee_inv, m_foo_inv;

  std::string m_mode;  //!< mult mode

  AFIELD m_v1;
  AFIELD m_v2;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  std::vector<Channel*> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  AFopr_CloverTerm<AFIELD>  *m_fopr_ct;

  Timer m_timer;

public:
  //! constructor (to be discarded).
  AFopr_Clover_eo(){ init(); }

  //! constructor.
  AFopr_Clover_eo(const Parameters& params){
    init();
    set_parameters(params);
  }

  //! destructor.
  ~AFopr_Clover_eo(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, real_t cSW, std::vector<int> bc);

  //! setting gauge configuration.
  void set_config(Field* u);
    
  void set_config(unique_ptr<Field_G>& u)
  { set_config(u.get()); }
    
  //! setting gauge configuration with boundary condition.
  void set_boundary_config(AFIELD&, int mu);

  //! returns the pointer to gauge configuration.
  inline Field* get_conf() { return m_conf; };

  //! setting mult mode.
  void set_mode(std::string mode){ m_mode = mode; }

  //! returns the current mult mode.
  std::string get_mode() const { return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_gm4(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
            const std::string mode);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void gm5(AFIELD&, const AFIELD&);

  void aypx(real_t, AFIELD&, const AFIELD&);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void mult_test(AFIELD&, const AFIELD&);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const int ieo);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo_alt(AFIELD&, const AFIELD&, const int ieo);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst2; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

 private:
  //! initial setup.
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  void mult_xp(real_t*, real_t*, const int);
  void mult_xm(real_t*, real_t*, const int);
  void mult_yp(real_t*, real_t*, const int);
  void mult_ym(real_t*, real_t*, const int);
  void mult_zp(real_t*, real_t*, const int);
  void mult_zm(real_t*, real_t*, const int);
  void mult_tp(real_t*, real_t*, const int);
  void mult_tm(real_t*, real_t*, const int);

  //! setting clover term inverse.
  void set_csw(int ieo);

  void mult_cswinv(AFIELD&, const AFIELD&, int ieo);

  void mult_cswinv(real_t *v2, real_t *v1, int ieo);

  void mult_cswinv_dirac(real_t *v2, real_t *v1,
                         real_t *csw_inv, int site);

  void mult_gm5(real_t*, real_t*);

  void mult_gm4(real_t*, real_t*);

  void clear(real_t*);
  void aypx(real_t, real_t*, real_t*);
  void scal(real_t*, const real_t);

  void clear(real_t*, const int, const int);
  void aypx(real_t, real_t*, real_t*, const int, const int);
  void scal(real_t*, const real_t, const int, const int);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Clover_eo",
                                               create_object_with_params);
    return init1;
  }
#endif

};

#endif  // AFOPR_WILSON_H
