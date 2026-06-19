/*!
        @file    corr2pt_Staggered.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef CORR2PT_STAGGERED_INCLUDED
#define CORR2PT_STAGGERED_INCLUDED

#include "Field/field_F_1spinor.h"
#include "Field/index_lex.h"
#include "Parameters/parameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Two-point correlator for staggered fermions.

/*!
   This is just a temporary implementation for test of
   staggered solver.
                              [28 Dec 2011 H.Matsufuru]
 */

class Corr2pt_Staggered
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Index_lex m_index;

 public:
  Corr2pt_Staggered()
    : m_vl(CommonParameters::Vlevel()) {}

  // optional
  Corr2pt_Staggered(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

 private:
  // non-copyable
  Corr2pt_Staggered(const Corr2pt_Staggered&);
  Corr2pt_Staggered& operator=(const Corr2pt_Staggered&);

 public:
  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void meson(
    std::vector<double>&,
    const std::vector<Field_F_1spinor>& sq1,
    const std::vector<Field_F_1spinor>& sq2);
};
#endif
