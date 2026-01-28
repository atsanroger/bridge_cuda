/*!
      @file    afopr_Domainwall_5din.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

#ifndef QXS_AFOPR_DOMAINWALL_5DIN_INCLUDED
#define QXS_AFOPR_DOMAINWALL_5DIN_INCLUDED

#include <vector>
#include <string>

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"


class Field;
class Field_G;

//! Optimal Domain-wall fermion operator.

/*!
   template version.                    [18 Apr 2017 H.Matsufuru]
 */
template<typename AFIELD>
class AFopr_Domainwall_5din : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                 //!< quark mass
  real_t m_M0;                 //!< domain-wall height
  int m_Ns;                    //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<real_t> m_b;     //!< coefficient b
  std::vector<real_t> m_c;     //!< coefficient c
  std::string m_repr;          //!< gamma matrix representation
  std::string m_kernel_type;   //!< kernel type (only Wilson is accepted)

  std::string m_mode;

  Bridge::VerboseLevel m_vl;  //!< verbose level

  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol, m_Ndim;
  int m_NinF, m_Nvcd, m_Ndf;
  int m_Nxv, m_Nyv, m_Nstv;

  AFIELD m_U;

  AFIELD m_w4, m_v4, m_t4, m_y4;   //!< working 4d vectors.
  AFIELD m_w1, m_v1, m_v2;         //!< working 5d vectors.

  // for preconditioning
  std::vector<real_t> m_dp;
  std::vector<real_t> m_dpinv;
  std::vector<real_t> m_dm;
  std::vector<real_t> m_e;
  std::vector<real_t> m_f;
  real_t m_g;

  int m_Nsize[4];  // lattice sizes in units of SIMD variable

  int do_comm[4];  //!< communication switch (4=Ndim): (0: n, 1: y).
  int do_comm_any; //!< communication switch (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  AFopr_Domainwall_5din(const Parameters& params)
    : AFopr<AFIELD>()
  {
    m_vl = CommonParameters::Vlevel();
    init(params);
  }

  ~AFopr_Domainwall_5din() { tidyup(); }

  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const double mq, const double M0,
                      const int Ns, const std::vector<int> bc,
                      const double b, const double c);

  void get_parameters(Parameters& params) const;

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<double> b,
                        const std::vector<double> c);

  //! this class needs convert of fermion field.
  bool needs_convert() { return true; }

  //! convert Field to AField for this class.
  void convert(AFIELD&, const Field&);

  //! reverse AField to Field.
  void reverse(Field&, const AFIELD&);

  //! setting gauge configuration (common interface).
  void set_config(Field *U);

  //  void set_boundary_config(AFIELD& U, const int mu);

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(AFIELD& v, const AFIELD& w)
  { mult(v, w, m_mode); }

  void mult_dag(AFIELD& v, const AFIELD& w)
  { mult_dag(v, w, m_mode); }

  //! mult with specified mode.
  //! Possible modes: D, Ddag, DdagD, D_prec, Ddag_prec, DdagD_prec.
  void mult(AFIELD& v, const AFIELD& w, std::string mode);

  //! mult_dag with specified mode.
  void mult_dag(AFIELD& v, const AFIELD& w, std::string mode);

  //! mult_dag with specified mode.
  void mult_gm5(AFIELD& v, const AFIELD& w);

  void DdagD(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);

  void mult_D(AFIELD&, const AFIELD&);
  void mult_Ddag(AFIELD&, const AFIELD&);

  //  preconditioner
  void DdagD_prec(AFIELD&, const AFIELD&);
  void D_prec(AFIELD&, const AFIELD&);
  void Ddag_prec(AFIELD&, const AFIELD&);
  void Prec(AFIELD&, const AFIELD&);
  void Precdag(AFIELD&, const AFIELD&);

  void L_inv(AFIELD&, const AFIELD&);
  void U_inv(AFIELD&, const AFIELD&);
  void Ldag_inv(AFIELD&, const AFIELD&);
  void Udag_inv(AFIELD&, const AFIELD&);

  int field_nin() { return m_NinF; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return 1; }

  //! this returns the number of floating point number operations.
  double flop_count() { return flop_count(m_mode); }

  //! flop-count for specified mode.
  double flop_count(std::string mode);

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! final tidyup.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! setting gauge configuration (setting omp parallel).
  void set_config_omp(Field *u);

  //! setting gauge configuration (implementation).
  void set_config_impl(Field *u);

  //! hopping part of fermion operator.
  void Dhop(real_t *, real_t *);

  void Dhop_1(real_t *, real_t *);

  void Dhop_2(real_t *, real_t *);

  void Dhop_b(real_t *, real_t *);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_5din(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
      "Domainwall_5din", create_object_with_params);
    init1 &= AFopr<AFIELD>::Factory_params::Register(
      "Domainwall_General_5din", create_object_with_params);
    return init1;
  }
#endif
};

// for transition of the class name
template<class AFIELD>
using AFopr_Domainwall_General_5din = AFopr_Domainwall_5din<AFIELD>;

#endif
