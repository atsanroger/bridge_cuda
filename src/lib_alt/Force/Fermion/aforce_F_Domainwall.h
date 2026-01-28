/*!
        @file    aforce_F_Domainwall.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-05-28 16:32:14 #$

        @version $LastChangedRevision: 2520 $
*/

#ifndef AFORCE_F_DOMAINWALL_INCLUDED
#define AFORCE_F_DOMAINWALL_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Fopr/afopr.h"
#include "lib/Fopr/afopr_Domainwall.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"


//! Force calculation for domain-wall fermions.

/*!
    The template version of domain-wall fermion operator.
                                [28 Dec 2011 H.Matsufuru]
 */


template<typename AFIELD>
class AForce_F_Domainwall : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                  //!< quark mass
  real_t m_M0;                  //!< domain-wall height
  int m_Ns;                     //!< size of fifth-dimension
  std::vector<int> m_boundary;  //!< boundary condition
  std::vector<real_t> m_b;      //!< coefficient b (array)
  std::vector<real_t> m_c;      //!< coefficient c (array)
  std::string m_repr;           //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;

  using AForce_F<AFIELD>::m_Ucp;

  std::string m_kernel_type;

  AFopr_Domainwall<AFIELD>* m_fopr_dw;
  AFopr<AFIELD>* m_fopr_w;
  AForce_F<AFIELD>* m_force_w;

  AFIELD m_v1, m_v2, m_v3;    //!< 4D working vector
  AFIELD m_z1;                //!< 5D working vector
  AFIELD m_force1, m_force2;  //!< working gauge field

 public:

  AForce_F_Domainwall(const Parameters& params)
    : AForce_F<AFIELD>()
  { init(params); }

  ~AForce_F_Domainwall() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const real_t mq,
                      const real_t M0,
                      const int Ns,
                      const std::vector<int> bc,
                      const real_t b,
                      const real_t c);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  void init(const Parameters& params);

  void tidyup();

  void set_config_omp(Field *U);

  void set_config_impl(Field *U);

};

#endif
