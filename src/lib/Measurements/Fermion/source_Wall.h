/*!
        @file    source_Wall.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SOURCE_WALL_INCLUDED
#define SOURCE_WALL_INCLUDED

#include "source.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wall source for 4-spinor fermion.

/*!
    This class sets a wall source vector
    for the 4-spinor (Wilson-type) fermion.
                                [02 Feb 2013 Y.Namekawa]
    Add set_all_color_spin,etc  [ 4 Apr 2017 Y.Namekawa]
 */

class Source_Wall : public Source
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Index_lex m_index;
  std::vector<int> m_source_position;
  bool m_in_node;

 public:
  Source_Wall() : m_vl(CommonParameters::Vlevel()) {}

  Source_Wall(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void get_parameters(Parameters& params) const;

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Wall();
  }

  static Source *create_object_with_params(const Parameters& params)
  {
    return new Source_Wall(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Source::Factory::Register("Wall", create_object);
    init &= Source::Factory_params::Register("Wall", create_object_with_params);
    return init;
  }
#endif
};
#endif /* SOURCE_WALL_INCLUDED */
