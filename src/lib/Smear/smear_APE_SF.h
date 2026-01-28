/*!
        @file    smear_APE_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SMEAR_APE_SF_INCLUDED
#define SMEAR_APE_SF_INCLUDED

#include "smear.h"

#include "lib/Measurements/Gauge/staple_SF.h"
#include "lib/Field/field_G.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! APE type smearing of link variables.

/*!
    This class is alternative to the Smear_APE class.
                            [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.    [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class Smear_APE_SF : public Smear
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim;
  std::valarray<double> m_rho;
  Projection *m_proj;

  std::vector<double> m_phi;   //!< SF boundary condition at t=0
  std::vector<double> m_phipr; //!< SF boundary condition at t=Nt

  Mat_SU_N m_wk;               //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;             //!< SF boundary condition at t=Nt

 public:
  Smear_APE_SF(Projection *proj)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc()),
    m_proj(proj) {}

  Smear_APE_SF(Projection *proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc()),
    m_proj(proj)
  { set_parameters(params); }

  ~Smear_APE_SF() {}

  void set_parameters(const Parameters& params);

  void set_parameters(const double rho1,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void set_parameters(const std::vector<double>& rho,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void smear(Field_G& Usmear, const Field_G& U);

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_APE_SF(proj);
  }

  static Smear *create_object_with_params(Projection *proj, const Parameters& params)
  {
    return new Smear_APE_SF(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory::Register("APE_SF", create_object);
    init &= Smear::Factory_params::Register("APE_SF", create_object_with_params);
    return init;
  }
#endif
};

#endif
