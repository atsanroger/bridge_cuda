/*!
        @file    source_Exponential.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SOURCE_EXPONENTIAL_INCLUDED
#define SOURCE_EXPONENTIAL_INCLUDED

#include "source.h"

#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Exponentially smeared source for 4-spinor fermion.

/*!
    This class sets an exponentially smeared source vector
    for the 4-spinor (Wilson-type) fermion.
                                      [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    Add set_all_color_spin,etc        [ 4 Apr 2017 Y.Namekawa]
 */

class Source_Exponential : public Source
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Index_lex m_index;
  std::vector<int> m_source_position;
  double m_slope, m_power;
  bool m_in_node;
  Field m_src_func;

 public:
  Source_Exponential() : m_vl(CommonParameters::Vlevel()) {}

  Source_Exponential(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position,
                      const double slope, const double power);

  void get_parameters(Parameters& params) const;

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Exponential();
  }

  static Source *create_object_with_params(const Parameters& params)
  {
    return new Source_Exponential(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Source::Factory::Register("Exponential", create_object);
    init &= Source::Factory_params::Register("Exponential", create_object_with_params);
    return init;
  }
#endif
};
#endif /* SOURCE_EXPONENTIAL_INCLUDED */
