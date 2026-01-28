/*!
        @file    aforceSmear_APE.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-03 23:32:09 #$

        @version $LastChangedRevision: 2529 $
*/

#ifndef QXS_AFORCESMEAR_APE_INCLUDED
#define QXS_AFORCESMEAR_APE_INCLUDED

#include "lib/Smear/aforceSmear.h"

#include "lib_alt_QXS/Field/shiftAField_lex.h"
#include "lib/Field/shiftField_lex.h"

#include "lib/Smear/aprojection.h"
#include "lib/Smear/smear_APE.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
                                [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

template<typename AFIELD>
class AForceSmear_APE : public AForceSmear<AFIELD>
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim, m_Nvol;
  std::vector<double> m_rho;
  AProjection<AFIELD>* m_proj;
  ShiftField_lex* m_shift;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  Field_G m_vt1, m_vt2, m_vt3;
  Field_G m_C, m_Ctmp, m_sigmap_tmp, m_Xi, m_sigma_tmp;

 public:
  AForceSmear_APE(AProjection<AFIELD> *proj)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj)
  {
    init();
  }

  AForceSmear_APE(AProjection<AFIELD>* proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj)
  {
    init();
    set_parameters(params);
  }

  ~AForceSmear_APE(){ tidyup(); }

  // Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  // Setting parameters with uniform smearing parameter.
  void set_parameters(const double rho1);

  // Setting parameters with anisotropic smearing parameter.
  void set_parameters(const std::vector<double>& rho);

  // Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  // Force computation.
  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

 private:
  void init();

  void tidyup();

  double rho(const int mu, const int nu)
  {
    return m_rho[mu + nu * m_Ndim];
  }

  void force_each(Field_G&,
                  const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);

  void staple(Field_G&,
              const Field_G&, const Field_G&,
              const int mu, const int nu);

#ifdef USE_FACTORY
 private:
  static AForceSmear<AFIELD> *create_object(AProjection<AFIELD> *proj)
  {
    return new AForceSmear_APE<AFIELD>(proj);
  }

  static AForceSmear<AFIELD> *create_object_with_params(
                     AProjection<AFIELD>* proj, const Parameters& params)
  {
    return new AForceSmear_APE<AFIELD>(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AForceSmear<AFIELD>::Factory::Register("APE", create_object);
    init &= AForceSmear<AFIELD>::Factory_params::Register("APE",
                                              create_object_with_params);
    return init;
  }
#endif
};
#endif
