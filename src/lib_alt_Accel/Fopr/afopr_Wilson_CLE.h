/*!
      @file    afopr_Wilson_CLE.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef ACCEL_AFOPR_WILSON_CLE_H
#define ACCEL_AFOPR_WILSON_CLE_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

class Field;

//! Wilson fermion operator for Complex Langevin equaiton.
/*!
    This class is a complexified standard Wilson fermion operator
    making use of the kernel code.
    CAUTION: since the link variables are out of SU(3), the
      reconstruction of the third column of matrix is not
      available. Remember to call only the kernel that directly
      load the original gauge field.
                                        [15 Jan 2021 H.Matsufuru]
*/

template<typename AFIELD>
class AFopr_Wilson_CLE : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  enum Gamma_repr {DIRAC, CHIRAL} m_repr; //!< gamma matrix representation

  int m_Nc, m_Nd, m_Nst, m_Ndim;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_chpot;             //!< chemical potential.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  Bridge::VerboseLevel m_vl;

  Field  *m_conf;       //!< poniter to the original gauge config.
  AFIELD  m_U;          //!< gauge config. (no boundary condition)
  AFIELD  m_Uinvd;      //!< (m_U^inv)^dag

  std::string  m_mode;  //!< mult mode

  AFIELD m_v2;     //!< working vector for DdagD
  AFIELD m_w1;     //!< working vector for Ddag

  int do_comm[4];  //!< switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; //!< switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  int m_Nsize[4];  //!< lattice sizes used in kernel code
  int m_bc[4];     //!< boundary condition used in kernel code
  int m_bc2[4];    //!< boundary condition in bulk used in kernel code

public:

  //! constructor.
  AFopr_Wilson_CLE(const Parameters& params) : AFopr<AFIELD>()
  { init();
    set_parameters(params); }

  //! destructor.
  ~AFopr_Wilson_CLE(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const real_t CKs,
                      const std::vector<int> bc, 
                      const real_t chpot);

  //! setting gauge configuration.
  void set_config(Field* u);

  //! setting gauge configuration (unique pointer).
  void set_config(unique_ptr<Field_G>& u){ set_config(u.get()); }

  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const { return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  //! mult with specified mode
  void mult(AFIELD&, const AFIELD&, const std::string mode);
  void mult_dag(AFIELD&, const AFIELD&, const std::string mode);

  //! upward nearest neighbor hopping term.
  virtual void mult_up(int mu, AFIELD&, const AFIELD&);

  //! downward nearest neighbor hopping term.
  virtual void mult_dn(int mu, AFIELD&, const AFIELD&);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

  //! returns floating operation counts.
  double flop_count(const std::string mode);

 private:

  //! initial setup.
  void init();

  //! final tidy-up.
  void tidyup();

  //! determine (U^\dag)^{-1}
  void inverse_dag(AFIELD&, AFIELD&);

  void mult_DdagD(AFIELD&, const AFIELD&);

  void mult_D(AFIELD&, const AFIELD&);
  void mult_D_alt(AFIELD&, const AFIELD&);
  void mult_Ddag(AFIELD&, const AFIELD&);
  void mult_Ddag_alt(AFIELD&, const AFIELD&);

  void mult_H(AFIELD&, const AFIELD&);

  //! setup channels for communication.
  void setup_channels();

  void mult_xp(real_t*, real_t*, int jd);
  void mult_xm(real_t*, real_t*, int jd);
  void mult_yp(real_t*, real_t*, int jd);
  void mult_ym(real_t*, real_t*, int jd);
  void mult_zp(real_t*, real_t*, int jd);
  void mult_zm(real_t*, real_t*, int jd);
  void mult_tp(real_t*, real_t*, int jd);
  void mult_tm(real_t*, real_t*, int jd);

  void mult_gm5(real_t*, real_t*);

  void clear(real_t*);
  void aypx(real_t, real_t*, real_t*);
  void axpy(real_t*, real_t, real_t*);
  void scal(real_t*, real_t);
  void gm5_aypx(real_t, real_t*, real_t*);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Wilson_CLE(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Wilson_CLE",
                                               create_object_with_params);
    return init1;
  }
#endif

};

#endif  // AFOPR_WILSON_CLE_H
