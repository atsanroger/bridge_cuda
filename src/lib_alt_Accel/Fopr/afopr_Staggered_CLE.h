/*!
      @file    afopr_Staggered_CLE.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef ACCEL_AFOPR_STAGGERED_CLE_INCLUDED
#define ACCEL_AFOPR_STAGGERED_CLE_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"


//! Staggered fermion operator for Complex Langevin equaiton.

/*!
    This class is a complexified standard staggered fermion operator.
    Note that this implementation adopt the mass normalization,
    not the hopping parameter normalization.
                                        [05 Mar 2020 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Staggered_CLE : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nc, m_Nvc, m_Nst, m_Ndim;

  real_t  m_mq;       //!< quark mass.
  real_t  m_chpot;    //!< chemical potential.
  std::vector<int>  m_boundary;  //!< boundary conditions.
  Bridge::VerboseLevel m_vl;     //!< verbose level

  AFIELD  m_staggered_phase;
  AFIELD  m_parity;
  AFIELD  m_U;      //!< gauge field multiplied by staggered phase.
  AFIELD  m_Uinvd;  //!< ((m_U^inv)^dag) multiplied by staggered phase.

  std::string m_mode;

  AFIELD m_v2;
  AFIELD m_w1;  //! working vector

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  int m_Nsize[4];
  int m_bc[4];
  int m_bc2[4];

 public:
  //! constructor.
  AFopr_Staggered_CLE(const Parameters& params) : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Staggered_CLE()
  { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const real_t mq,
                      const std::vector<int> bc,
                      const real_t chpot);

  //! setting gauge configuration.
  void set_config(Field *U);

  void set_config(unique_ptr<Field_G>& U)
  { set_config(U.get()); }

  void set_mode(std::string mode){ m_mode = mode; }

  std::string get_mode() const {  return m_mode; }


  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);


  void mult_gm5(AFIELD&);

  void mult_staggered_phase(AFIELD&, int mu);

  void fprop_normalize(AFIELD& v) { return; }

  void fopr_normalize(AFIELD& v) { return; }

  int field_nvol() { return CommonParameters::Nvol(); }
  int field_nin() { return 2 * CommonParameters::Nc(); }
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

  //! returns floating operation counts.
  double flop_count(const std::string mode);

 private:
  void init(const Parameters& params);
  void tidyup();

  void set_staggered_phase();
  //  void set_staggered_config(AFIELD_G *Ueo);

  void setup_channels();

  void inverse_dag(AFIELD&, AFIELD&);

  void D(AFIELD&, const AFIELD&);
  void D_alt(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void Ddag_alt(AFIELD&, const AFIELD&);
  void DdagD(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);

  void clear(real_t*);
  void aypx(real_t, real_t*, real_t*);
  void axpy(real_t*, real_t, real_t*);
  void scal(real_t*, real_t);

  void mult_xp(real_t*, real_t*, int jd);
  void mult_xm(real_t*, real_t*, int jd);
  void mult_yp(real_t*, real_t*, int jd);
  void mult_ym(real_t*, real_t*, int jd);
  void mult_zp(real_t*, real_t*, int jd);
  void mult_zm(real_t*, real_t*, int jd);
  void mult_tp(real_t*, real_t*, int jd);
  void mult_tm(real_t*, real_t*, int jd);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Staggered_CLE(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Staggered_CLE",
                                               create_object_with_params);
    return init1;
  }
#endif

};

#endif
