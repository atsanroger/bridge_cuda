/*!
        @file    source_Staggered_Wall_Cube.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOURCE_STAGGERED_WALL_CUBE_INCLUDED
#define SOURCE_STAGGERED_WALL_CUBE_INCLUDED

#include "lib/Measurements/Fermion/source.h"
#include "lib/Parameters/parameters.h"
#include "lib/Field/index_lex.h"
#include "lib/Field/field_F_1spinor.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;


//! Setting source vector with wall (cube) source for staggered fermion.

/*!
    This class set the wall source with one of corner site of 2x2x2
    spatial cube is set to one.
    The nonzero site (ix,iy,iz) of the cube is specified by isrc, as
    isrc = ix + 2*iy + 4*iz and this in [0,7], given in set() method.
    The value of isrc may be provide through idx = ic + Nc * isrc.
    This cubic source is repeatedly applied to the space at time t
    specified by parameter source_position.
    Gauge fixing is needed before the measurement.
                                             [15 Oct 2020 H.Matsufuru]
    Modified so that odd values of Nx, Ny, Nz are allowed.
                                             [14 Jan 2023 H.Matsufuru]
 */

class Source_Staggered_Wall_Cube : public Source
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_t_src;         //!< source time slice.
  Index_lex m_index;

 public:
  //! constructor
  Source_Staggered_Wall_Cube() : Source() { init(); }

 private:
  // non-copyable
  Source_Staggered_Wall_Cube(const Source_Staggered_Wall_Cube&);
  Source_Staggered_Wall_Cube& operator=(const Source_Staggered_Wall_Cube&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const int source_position);

  void get_parameters(Parameters& params) const;

  void set(Field& src, const int idx);
  void set(Field& src, const int ic, const int isrc);

  void set(Field_F_1spinor& src, const int ic, const int i_src);

  void set_all_color(Field&, int)
  {
    vout.crucial("%s: set_all_color() is not available.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  void set_all_color_spin(Field&)
  {
    vout.crucial("%s: set_all_color_spin() is not available.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

 private:
  void init();


#ifdef USE_FACTORY
 private:
  static Source *create_object()
  { return new Source_Staggered_Wall_Cube(); }

 public:
  static bool register_factory()
  {
    return Source::Factory::Register("Staggered_Wall_Cube", create_object);
  }
#endif
};
#endif
