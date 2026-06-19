/*!
      @file    afopr_Staggered_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2554 $
*/

#ifndef ACCEL_AFOPR_STAGGERED_EO_INCLUDED
#define ACCEL_AFOPR_STAGGERED_EO_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;


//! Staggered fermion operator.

/*!
    This class is a complexified standard staggered fermion operator.
    The implementation is based on the even-odd staggered version.
                                        [19 Sep 2015 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Staggered_eo : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nc, m_Nvc, m_Nst, m_Ndim;
  int m_Nst2;
  int m_Ieo_origin;

  real_t  m_mq;                  //!< quark mass.
  std::vector<int>  m_boundary;  //!< boundary conditions.
  Bridge::VerboseLevel m_vl;     //!< verbose level

  AFIELD  m_staggered_phase;
  AFIELD  m_parity;
  AFIELD  m_Ueo;    //!< gauge field multiplied by staggered phase.

  std::string m_mode;

  AFIELD m_v1;
  AFIELD m_v2;

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
  AFopr_Staggered_eo(const Parameters& params) : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Staggered_eo()
  { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const real_t mq, const std::vector<int> bc);

  //! getting parameters via a Parameters object.
  void get_parameters(Parameters& params) const;

  //! setting gauge configuration.
  void set_config(Field *U);

  void set_mode(std::string mode);

  std::string get_mode() const {  return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&, const std::string mode);

  void mult_gm5(AFIELD&);

  void mult_staggered_phase(AFIELD&, int mu);

  void normalize_fprop(AFIELD& v);

  void normalize_fopr(AFIELD& v);

  int field_nvol(){ return m_Nst2; }
  int field_nin(){  return m_Nvc; }
  int field_nex(){  return 1; }

  //! returns floating operation counts.
  double flop_count();

  //! returns floating operation counts.
  double flop_count(const std::string mode);

 private:
  void init(const Parameters& params);
  void tidyup();

  void set_staggered_phase();

  void setup_channels();

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const int ieo);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const AFIELD&,
           const int ieo, const int iflag);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo_alt(AFIELD&, const AFIELD&, const int ieo);

  void D_alt(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void DdagD(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);

  void scal(AFIELD&, real_t);

  void clear(real_t*);
  void aypx(real_t, real_t*, real_t*);
  void scal(real_t*, real_t);

  void mult_xp(real_t*, real_t*, const int);
  void mult_xm(real_t*, real_t*, const int);
  void mult_yp(real_t*, real_t*, const int);
  void mult_ym(real_t*, real_t*, const int);
  void mult_zp(real_t*, real_t*, const int);
  void mult_zm(real_t*, real_t*, const int);
  void mult_tp(real_t*, real_t*, const int);
  void mult_tm(real_t*, real_t*, const int);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Staggered_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Staggered_eo",
                                                 create_object_with_params);
    return init1;
  }
#endif

};

#endif
