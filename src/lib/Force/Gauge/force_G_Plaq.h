/*!
        @file    force_G_Plaq.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCE_G_PLAQ_INCLUDED
#define FORCE_G_PLAQ_INCLUDED

#include "force_G.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC force class for plaquette gauge action.

/*!
    Standard plaquette gauge action.
                             [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
 */

class Force_G_Plaq : public Force_G
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  //- NB. m_U has been defined in force_G.h
  // Field_G *m_U;

  double m_beta;
  Staple_lex m_staple;

 public:
  Force_G_Plaq() : m_vl(CommonParameters::Vlevel()) {}

  Force_G_Plaq(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  ~Force_G_Plaq() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta);

  void get_parameters(Parameters& params) const;

  //- NB. set_config has been defined in force_G.h
  // void set_config(Field *U);

  void force_core(Field& force);

#ifdef USE_FACTORY
 private:
  static Force_G *create_object()
  {
    return new Force_G_Plaq();
  }

  static Force_G *create_object_with_params(const Parameters& params)
  {
    return new Force_G_Plaq(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Force_G::Factory::Register("Force_G_Plaq", create_object);
    init &= Force_G::Factory_params::Register("Force_G_Plaq", create_object_with_params);
    return init;
  }
#endif
};
#endif
