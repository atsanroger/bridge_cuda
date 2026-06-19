/*!
        @file    smear_APE.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SMEAR_APE_INCLUDED
#define SMEAR_APE_INCLUDED

#include "smear.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! APE type smearing of link variables.

/*!
    This class smears link variables with APE-type construction
    of smeared links with a given projection operator to SU(N)
    group element.
    Parameter is \rho(\mu,\nu), which in general depends on
    the directions of the modified link and the staple.
    By explicitly giving \rho(\mu,\nu) as std::vector object,
    anisotropic setup is possible, while isotropic setup
    requires only one double parameter, `rho_uniform'.
                       [08 Apr 2012/15 Jul 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

class Smear_APE : public Smear
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim;                    //!< spacetime dimension
  std::valarray<double> m_rho;   //!< smearing parameter
  Projection *m_proj;            //!< projector to group element.

 public:
  //! Constructor requires a pointer to Projection object.
  Smear_APE(Projection *proj)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_proj(proj) {}

  Smear_APE(Projection *proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_proj(proj)
  {
    set_parameters(params);
  }

  //! Deconstructor
  ~Smear_APE() {}

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting parameter with isotropic parameter.
  void set_parameters(const double rho1);

  //! Setting parameter with anisotropic parameter.
  void set_parameters(const std::vector<double>& rho);

  //! Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

 private:

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_APE(proj);
  }

  static Smear *create_object_with_params(Projection *proj, const Parameters& params)
  {
    return new Smear_APE(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory::Register("APE", create_object);
    init &= Smear::Factory_params::Register("APE", create_object_with_params);
    return init;
  }
#endif
};
#endif
