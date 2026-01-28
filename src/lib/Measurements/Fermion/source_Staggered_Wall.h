/*!
        @file    source_Staggered_Wall.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SOURCE_STAGGERED_WALL_INCLUDED
#define SOURCE_STAGGERED_WALL_INCLUDED

#include "Field/field_F_1spinor.h"
#include "Field/index_lex.h"
#include "Parameters/parameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Setting source vector with wall source for staggered fermion.

/*!
    This is a temporary implementation for test of staggered
    fermion.
                                 [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
 */

class Source_Staggered_Wall
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_t_src;
  Index_lex m_index;

  Field src_wall_1;
  Field src_wall_2;

 public:
  Source_Staggered_Wall()
    : m_vl(CommonParameters::Vlevel()) {}

  Source_Staggered_Wall(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

 private:
  // non-copyable
  Source_Staggered_Wall(const Source_Staggered_Wall&);
  Source_Staggered_Wall& operator=(const Source_Staggered_Wall&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const int source_position);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void get_parameters(Parameters& params) const;

  void set(Field_F_1spinor& src, const int ic, const int i_src);

 private:
  void init();
};
#endif
