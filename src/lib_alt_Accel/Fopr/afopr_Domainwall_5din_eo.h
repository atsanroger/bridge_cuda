/*!
        @file    afopr_Domainwall_5din_eo.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-01-16 21:43:34 #$
        @version $LastChangedRevision: 2631 $
*/

#ifndef ACCEL_AFOPR_DOMAINWALL_5DIN_EO_INCLUDED
#define ACCEL_AFOPR_DOMAINWALL_5DIN_EO_INCLUDED

#include <vector>
#include <string>

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
#include "lib/Tools/timer.h"

using Bridge::vout;

#include "lib_alt_Accel/Field/aindex_eo.h"

class Field;

//! Domain-wall fermion operator with even-odd site index.

/*!
   This class implements the even-odd version of Domain-wall
   fermion including the Mobius form.
                                        [18 Apr 2017 H.Matsufuru]
   This class was developed as a template class in the alternative
   branch and then incorporated into ver.2.0.
                                        [07 Mar 2022 H.Matsufuru]
   Neff's improving parameter alpha is incorporated.
                                        [13 Jan 2025 H.Matsufuru]
 */
template<typename AFIELD>
class AFopr_Domainwall_5din_eo : public AFopr<AFIELD>
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
  std::vector<real_t> m_b;     //!< coefficient b (array)
  std::vector<real_t> m_c;     //!< coefficient c (array)
  real_t m_alpha;              //!< parameter alpha
  Bridge::VerboseLevel m_vl;   //!< verbose level
  std::string m_kernel_type;   //!< fermion kernel type
  std::string m_repr; 	       //!< Dirac matrix representation
  std::string m_impl; 	       //!< 4d or 5d kernel implementation

  std::vector<real_t> m_mat_inv;  //!< 5d matrix inverse

  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol, m_Ndim;
  int m_Nx2, m_Nst2;
  int m_NinF, m_Nvcd, m_Ndf, m_Nd, m_Nd2;
  int m_Ieo_origin;

  AIndex_eo<real_t,AFIELD::IMPL> m_index_eo;

  AFIELD m_Ueo;   //!< copied gauge config. with boundary conditions.

  AFIELD m_w1, m_v1, m_v2;        //!< woking vectors

  // for preconditioning
  std::vector<real_t> m_dp;
  std::vector<real_t> m_dpinv;
  std::vector<real_t> m_dm;
  std::vector<real_t> m_e;
  std::vector<real_t> m_f;
  real_t m_g;

  int do_comm[4];  //!< communication switch (4=Ndim): (0: n, 1: y)
  int do_comm_any; //!< communication switch (if any): (0: n, 1: y)

  int m_Nsize[4];  // lattice sizes in units of SIMD variable
  int m_bc[4];
  int m_bc2[4];

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  std::string m_mode;

 public:
  //! constructor.
  AFopr_Domainwall_5din_eo(const Parameters& params)
    : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Domainwall_5din_eo() { tidyup(); }

  void setup_channels();

  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const real_t b, const real_t c,
                      const real_t alpha);

  void get_parameters(Parameters& params) const;

  //! set parameters of kernel operaotr.
  //  void set_kernel_parameters(const Parameters& params);

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<real_t> b,
                        const std::vector<real_t> c);

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set parameters for preconditioning.
  void set_matrix5d_inverse();

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

  void mult(AFIELD& v, const AFIELD& w, const std::string mode);

  void mult_dag(AFIELD& v, const AFIELD& w, std::string mode);

  void mult_gm5(AFIELD& v, const AFIELD& w);

  void DdagD(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  // void H(AFIELD&, const AFIELD&);
  // void Hdag(AFIELD&, const AFIELD&);

  void DdagD_alt(AFIELD&, const AFIELD&);
  void D_alt(AFIELD&, const AFIELD&);
  void Ddag_alt(AFIELD&, const AFIELD&);

  void D_ee(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_ee(AFIELD&, const AFIELD&, const int ieo);

  void D_eo(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_eo(AFIELD&, const AFIELD&, const int ieo);

  void D_ee_inv(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_ee_inv(AFIELD&, const AFIELD&, const int ieo);

  void D_ee_inv_alt(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_ee_inv_alt(AFIELD&, const AFIELD&, const int ieo);

  void mult_gm5R(AFIELD&, const AFIELD&);
  void mult_R(AFIELD&, const AFIELD&);

  //  preconditioner
  void LU_inv(AFIELD&, const AFIELD&);
  void LUdag_inv(AFIELD&, const AFIELD&);

  int field_nin() { return m_NinF; }
  int field_nvol() { return m_Nst2; }
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

  int mat_index(int ida, int isa, int idb, int isb){
    return ida + m_Nd2 * isa + (m_Nd2 * m_Ns) * (idb + m_Nd2 * isb);
  }


  unique_ptr<Timer> timer_pack;
  unique_ptr<Timer> timer_bulk;
  unique_ptr<Timer> timer_boundary;
  unique_ptr<Timer> timer_comm;
  unique_ptr<Timer> timer_comm_recv_wait;
  unique_ptr<Timer> timer_comm_send_wait;
  unique_ptr<Timer> timer_comm_send_start;
  unique_ptr<Timer> timer_comm_recv_start;
  unique_ptr<Timer> timer_comm_test_all;
  unique_ptr<Timer> timer_mult_Deo;
  unique_ptr<Timer> timer_mult_Dee_inv;
  unique_ptr<Timer> timer_clear;

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_5din_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
      "Domainwall_5din_eo", create_object_with_params);
    return init1;
  }
#endif
};

#endif
