/*!
        @file    smear_APE_spatial.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SMEAR_APE_SPATIAL_INCLUDED
#define SMEAR_APE_SPATIAL_INCLUDED

#include "smear.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! APE type smearing of spatial link variables.

/*!
    This class smears spatial link variables with APE-type
    construction of smeared links with a given projection
    operator to SU(N) group element.
    Parameter is \rho, which specifies the mixing rate
    of original thin link and staples.
                            [09 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.    [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class Smear_APE_spatial : public Smear
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim;          //!< spacetime dimension
  double m_rho;        //!< smearing parameter
  Projection *m_proj;  //!< projector to group element.

 public:
  //! Constructor requires a pointer to Projection object.
  Smear_APE_spatial(Projection *proj)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()),
    m_proj(proj) {}

  Smear_APE_spatial(Projection *proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_Ndim(CommonParameters::Ndim()),
    m_proj(proj)
  {
    set_parameters(params);
  }

  //! Deconstructor
  ~Smear_APE_spatial() {}

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting smearing parameter.
  void set_parameters(const double rho);

  //! Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

 private:

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_APE_spatial(proj);
  }

  static Smear *create_object_with_params(Projection *proj, const Parameters& params)
  {
    return new Smear_APE_spatial(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory::Register("APE_spatial", create_object);
    init &= Smear::Factory_params::Register("APE_spatial", create_object_with_params);
    return init;
  }
#endif
};
#endif
