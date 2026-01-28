/*!
        @file    aprojection_Stout_SU3.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef APROJECTION_STOUT_SU3_INCLUDED
#define APROJECTION_STOUT_SU3_INCLUDED

#include "lib/Smear/aprojection.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Stout(exponential)-type projection to SU(N) gauge group.

/*!
    Present implementation applies to SU(3) case only.
    The SU(3) properties are explicitly used.
                                    [08 Apr 2012 H.Matsufuru]
    Add a strict check for Nc, complement to assert.
                                    [31 May 2021 Y.Namekawa]
 */

template<typename AFIELD>
class AProjection_Stout_SU3 : public AProjection<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  typedef typename AFIELD::complex_t complex_t;
  static const std::string class_name;

 private:
  int m_Ndf, m_Nst;
  int m_Nstv;
  Bridge::VerboseLevel m_vl;

  unsigned long int m_flop;
  double m_time;

  AFIELD m_v1, m_v2, m_v3;

 public:
  AProjection_Stout_SU3() { init(); }

  AProjection_Stout_SU3(const Parameters& params) {
    init();
    set_parameters(params);
  }

  ~AProjection_Stout_SU3() {}

  void set_parameters(const Parameters& param);

  void get_parameters(Parameters& param) const;

  //! projection U = P[alpha, C, Uorg]
  void project(Field_G& U,
               const double alpha,
               const Field_G& Cst, const Field_G& Uorg);

  //! temporary implementation
  void project_alt(Field_G& U,
               const double alpha,
               const Field_G& Cst, const Field_G& Uorg);

  //! projection with AFIELD
  void project(AFIELD& U, const double alpha,
               const AFIELD& Cst, const AFIELD& Uorg);



  //! determination of fields for force calculation
  void force_recursive(Field_G& Xi, Field_G& iTheta,
                       const double alpha, const Field_G& Sigmap,
                       const Field_G& C, const Field_G& U);

  void print_stat();

 private:
  void init();

  void exp_iQ(Field_G& e_iQ, const Field_G& iQ);
  void exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ);

  void set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
              const double& u, const double& w);

  void set_uw(double& u, double& w,
              const Mat_SU_N& iQ2, const Mat_SU_N& iQ3);

  double func_xi0(const double w);
  double func_xi1(const double w);

#ifdef USE_FACTORY
 private:
  static AProjection<AFIELD> *create_object()
  {
    return new AProjection_Stout_SU3<AFIELD>();
  }

  static AProjection<AFIELD> *create_object_with_params(
                                             const Parameters& params)
  {
    return new AProjection_Stout_SU3<AFIELD>(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AProjection<AFIELD>::Factory::Register("Stout_SU3", create_object);
    init &= AProjection<AFIELD>::Factory_params::Register("Stout_SU3",
                                            create_object_with_params);
    return init;
  }
#endif
};
#endif
