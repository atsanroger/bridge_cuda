/*!
        @file    source_Staggered_Wall_Evenodd.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOURCE_STAGGERED_WALL_EVENODD_INCLUDED
#define SOURCE_STAGGERED_WALL_EVENODD_INCLUDED

#include "source.h"
#include "Field/field_F_1spinor.h"
#include "Field/index_lex.h"
#include "Parameters/parameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Setting source vector with wall source for staggered fermion.

/*!
    This class set the wall source for the staggered fermion spectroscopy.
    On give time slice, the source vector is set to 1 for all spacial
    sites (even source) or (-1)^{x+y+z} (odd source).
    Gauge fixing is needed before the measurement.
    Cf. R.Gupta et al., Phys. Rev. D 43, 2003 (1991).
                          [28 Dec 2011, updated 15 Oct 2020 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                        [14 Nov 2012 Y.Namekawa]
    Renamed: Source_Staggered_Wall -> Source_Staggered_Wall_Evenodd.
    Implementation is also slightly modified.
                                               [12 Oct 2020 H.Matsufuru]
   Modified so that odd values of Nx, Ny, Nz are allowed.
                                               [14 Jan 2023 H.Matsufuru]
 */

class Source_Staggered_Wall_Evenodd : public Source
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_t_src;         //!< source time slice
  Index_lex m_index;

 public:

  Source_Staggered_Wall_Evenodd() : Source() { init(); }

 private:
  // non-copyable
  Source_Staggered_Wall_Evenodd(const Source_Staggered_Wall_Evenodd&);
  Source_Staggered_Wall_Evenodd& operator=(const Source_Staggered_Wall_Evenodd&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const int source_position);

  void get_parameters(Parameters& params) const;

  void set(Field& src, const int idx);
  void set(Field& src, const int ic, const int isrc);

  void set(Field_F_1spinor& src, const int ic, const int isrc);

  void set_all_color(Field&, int)
  {
    vout.crucial("%s: function set_all_color is unavailable.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  void set_all_color_spin(Field&)
  {
    vout.crucial("%s: function set_all_color_spin is unavailable.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

 private:

  void init();


#ifdef USE_FACTORY
 private:
  static Source *create_object()
  { return new Source_Staggered_Wall_Evenodd(); }

 public:
  static bool register_factory()
  {
    return Source::Factory::Register("Staggered_Wall_Evenodd", create_object);
  }
#endif
};
#endif
