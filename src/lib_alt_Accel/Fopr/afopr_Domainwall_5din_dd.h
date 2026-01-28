/*!
      @file    afopr_Domainwall_5din_dd.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2587 $
*/

#ifndef ACCEL_AFOPR_DOMAINWALL_5DIN_DD_INCLUDED
#define ACCEL_AFOPR_DOMAINWALL_5DIN_DD_INCLUDED

#include <vector>
#include <string>

//#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Fopr/afopr_dd.h"

class Field;
class Field_G;

//! Domain-wall fermion operator.

/*!
    This class implements the domain-wall fermion operator.
    The first version was implemented only for the standard
    (Shamir's) form.                [24 Dec 2011 H.Matsufuru]
    Later it was generalized to cover the Mobius form as a
    template class in the alternative branch.
                                    [26 Feb 2019 H.Matsufuru]
    In ver.2.0, it was merged to trunk by renaming to
    AFopr_Domainwall.
                                    [06 Mar 2022 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Domainwall_5din_dd : public AFopr_dd<AFIELD>
{
 public:
  typedef typename AFIELD::real_t  real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                 //!< quark mass
  real_t m_M0;                 //!< domain-wall height
  int m_Ns;                    //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<real_t> m_b;     //!< coefficient b (array)
  std::vector<real_t> m_c;     //!< coefficient c (array)
  Bridge::VerboseLevel m_vl;   //!< verbose level
  std::string m_kernel_type;
  std::string m_repr;

  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol, m_Ndim;
  int m_NinF, m_Nvcd, m_Ndf;

  // for preconditioning
  std::vector<real_t> m_dp;
  std::vector<real_t> m_dpinv;
  std::vector<real_t> m_dm;
  std::vector<real_t> m_e;
  std::vector<real_t> m_f;
  real_t m_g;

  int m_Nsize[4];  // lattice sizes in units of SIMD variable
  int m_bc[4];
  int m_bc2[4];

  int do_comm[4];  //!< communication switch (4=Ndim): (0: n, 1: y)
  int do_comm_any; //!< communication switch (if any): (0: n, 1: y)

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  std::vector<int> m_block_size;  //!< block size

  int m_Ieo;   //!< even-odd label of origin in units of block

  AFIELD m_Ublock;  //!< copied gauge config. with block condition.
  AFIELD m_Uhop;    //!< = m_U - m_Ublock

  std::string m_mode;

  AFIELD m_U;

  //  AFIELD m_w4, m_v4, m_t4, m_y4; //!< working 4d vectors.
  AFIELD m_w1, m_v1, m_v2;       //!< working 5d vectors.
  AFIELD m_w2, m_w3, m_w4;       //!< additional working vector.

  // for convert and reverse
  Field m_w4lex;
  AFIELD m_v4lex;


 public:
  //! constructor.
  AFopr_Domainwall_5din_dd(const Parameters& params) : AFopr_dd<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Domainwall_5din_dd() { tidyup(); }

  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const real_t b, const real_t c,
		      const std::vector<int> block_size);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const std::vector<real_t> vec_b,
                      const std::vector<real_t> vec_c,
		      const std::vector<int> block_size);

  void get_parameters(Parameters& params) const;

  void setup_channels();

  //! set parameters of kernel operaotr.
  //  void set_kernel_parameters(const Parameters& params);

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<real_t> b,
                        const std::vector<real_t> c);

  //! this class needs convert of fermion field.
  bool needs_convert() { return true; }

  //! convert Field to AField for this class.
  void convert(AFIELD&, const Field&);

  //! reverse AField to Field.
  void reverse(Field&, const AFIELD&);

  void set_config(Field *U);

  void set_config_omp(Field *U);

  void set_config_impl(Field *U);

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w);

  //! mult with specified mode.
  //! Possible modes are: "D", "Ddag", "DdagD", "DDdag",
  //!   "Prec", "D_prec", "Ddag_prec", "DdagD_prec".
  void mult(AFIELD& v, const AFIELD& w, std::string mode);

  //! mult_dag with specified mode.
  void mult_dag(AFIELD& v, const AFIELD& w, std::string mode);

  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_chproj_4d(AFIELD&, const AFIELD&, const int ipm);

  void mult_sap(AFIELD&, const AFIELD&, const int ieo);

  void mult_D_sap(AFIELD&, const AFIELD&, const int ieo);

  void mult_Ddag_sap(AFIELD&, const AFIELD&, const int ieo);

  void mult_dd(AFIELD&, const AFIELD&);

  void mult_D_dd(AFIELD&, const AFIELD&);

  void mult_Ddag_dd(AFIELD&, const AFIELD&);

  void mult_dup(AFIELD&, const AFIELD&, const int mu);

  void mult_ddn(AFIELD&, const AFIELD&, const int mu);


  int field_nin() { return m_NinF; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return 1; }

  //! this returns the number of floating point number operations.
  double flop_count() { return flop_count(m_mode); }

  //! flop-count for specified mode.
  double flop_count(std::string mode);

  //! flop-count for SAP.
  double flop_count_sap();

  void DdagD(AFIELD&, const AFIELD&);
  void DDdag(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void Hdag(AFIELD&, const AFIELD&);
  void mult_gm5R(AFIELD&, const AFIELD&);
  void mult_R(AFIELD&, const AFIELD&);

  //  preconditioner
  void DdagD_prec(AFIELD&, const AFIELD&);
  void D_prec(AFIELD&, const AFIELD&);
  void Ddag_prec(AFIELD&, const AFIELD&);
  void Prec(AFIELD&, const AFIELD&);
  void Precdag(AFIELD&, const AFIELD&);

  void D_alt(AFIELD&, const AFIELD&);
  void Ddag_alt(AFIELD&, const AFIELD&);

  void mult_up(AFIELD&, const AFIELD&, int isap, int idag, int mu);

  void mult_dn(AFIELD&, const AFIELD&, int isap, int idag, int mu);

 private:

  //! initial setup.
  void init(const Parameters& params);

  //! final tidyup.
  void tidyup();

  void L_inv(AFIELD&, const AFIELD&);
  void U_inv(AFIELD&, const AFIELD&);
  void Ldag_inv(AFIELD&, const AFIELD&);
  void Udag_inv(AFIELD&, const AFIELD&);

  void mult_gm5_4d(AFIELD&, const AFIELD&);

  void mult_up4(real_t*, real_t*, int isap, int idag, int mu);

  void mult_dn4(real_t*, real_t*, int isap, int idag, int mu);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_5din_dd(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AFopr<AFIELD>::Factory_params::Register(
                     "Domainwall_5din_dd", create_object_with_params);
    return init;
  }
#endif
};

#endif
