/*!
        @file    action_G_Rectangle.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef ACTION_G_RECTANGLE_INCLUDED
#define ACTION_G_RECTANGLE_INCLUDED

#include "Action/action.h"
#include "Force/Gauge/force_G_Rectangle.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC action class for rectangular gauge action.

/*!
    Gauge action with plaquette and rectangular Wilson loops.
    Iwasaki, Luscher-Weisz, DBW2 are examples of this type
    of action.
                                   [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
 */


class Action_G_Rectangle : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_beta;
  double m_c_plaq;
  double m_c_rect;
  std::string m_label;

  Field_G *m_U;
  Staple_lex m_staple;
  ShiftField_lex m_shift;
  Force_G *m_force_G;

 public:
  Action_G_Rectangle()
    : m_vl(CommonParameters::Vlevel())
  {
    m_force_G = Force_G::New("Force_G_Rectangle");
  }

  Action_G_Rectangle(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    // m_force_G = Force_G::New("Force_G_Rectangle");
    m_force_G = new Force_G_Rectangle(params);
    set_parameters(params);
  }

  ~Action_G_Rectangle()
  {
    delete m_force_G;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta,
                      const double c_plaq, const double c_rect);

  void get_parameters(Parameters& params) const;

  void set_label(const std::string label)
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
    return new Action_G_Rectangle();
  }

  static Action *create_object_with_params(const Parameters& params)
  {
    return new Action_G_Rectangle(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Action::Factory::Register("Action_G_Rectangle", create_object);
    init &= Action::Factory_params::Register("Action_G_Rectangle", create_object_with_params);
    return init;
  }
#endif
};
#endif
