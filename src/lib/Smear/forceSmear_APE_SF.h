/*!
        @file    forceSmear_APE_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCESMEAR_APE_ALT_INCLUDED
#define FORCESMEAR_APE_ALT_INCLUDED

#include "forceSmear.h"

#include "Field/field_SF.h"
#include "Field/shiftField_lex.h"

#include "Smear/smear_APE_SF.h"
#include "Smear/projection.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
                                [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */


class ForceSmear_APE_SF : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim, m_Nvol;
  std::vector<double> m_rho;
  Projection *m_proj;
  ShiftField_lex m_shift;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  //! SF boundary condition at t=0
  std::vector<double> m_phi;
  //! SF boundary condition at t=Nt
  std::vector<double> m_phipr;

  Mat_SU_N m_wk;    //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;  //!< SF boundary condition at t=Nt

  //  Field_G_SF set_wk;

 public:
  ForceSmear_APE_SF(Projection *proj)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  {
    init();
  }

  ForceSmear_APE_SF(Projection *proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  {
    init();
    set_parameters(params);
  }

  //  ~ForceSmear_APE_SF(){
  //  };

  void set_parameters(const Parameters& params);

  void set_parameters(const double rho1,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void set_parameters(const std::vector<double>& rho,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

  // old interface
  //Field force_udiv(const Field_G& Sigma, const Field_G& U);

 private:
  void init();

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
  static ForceSmear *create_object(Projection *proj)
  {
    return new ForceSmear_APE_SF(proj);
  }

  static ForceSmear *create_object_with_params(Projection *proj, const Parameters& params)
  {
    return new ForceSmear_APE_SF(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= ForceSmear::Factory::Register("APE_SF", create_object);
    init &= ForceSmear::Factory_params::Register("APE_SF", create_object_with_params);
    return init;
  }
#endif
};
#endif
