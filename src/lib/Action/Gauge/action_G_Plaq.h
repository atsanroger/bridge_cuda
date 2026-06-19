/*!
        @file    action_G_Plaq.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef ACTION_G_PLAQ_INCLUDED
#define ACTION_G_PLAQ_INCLUDED

#include "Action/action.h"
#include "Force/Gauge/force_G_Plaq.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action class for plaquette gauge action.

/*!
    Standard plaquette gauge action.
                             [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
 */


class Action_G_Plaq : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_beta;
  std::string m_label;

  Field_G *m_U;
  Staple_lex m_staple;
  Force_G *m_force_G;

 public:
  Action_G_Plaq()
    : m_vl(CommonParameters::Vlevel())
  {
    m_force_G = Force_G::New("Force_G_Plaq");
  }

  Action_G_Plaq(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    // m_force_G = Force_G::New("Force_G_Plaq");
    m_force_G = new Force_G_Plaq(params);
    set_parameters(params);
  }

  ~Action_G_Plaq()
  {
    delete m_force_G;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta);

  void get_parameters(Parameters& params) const;

  void set_label(std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  {
    return m_label;
  }

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

#ifdef USE_FACTORY
 private:
  static Action *create_object()
  {
    return new Action_G_Plaq();
  }

  static Action *create_object_with_params(const Parameters& params)
  {
    return new Action_G_Plaq(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Action::Factory::Register("Action_G_Plaq", create_object);
    init &= Action::Factory_params::Register("Action_G_Plaq", create_object_with_params);
    return init;
  }
#endif
};
#endif
