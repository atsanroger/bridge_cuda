/*!
        @file    fprop_alt.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-11-13 13:02:59 #$
        @version $LastChangedRevision: 2553 $
*/

#ifndef FPROP_ALT_INCLUDED
#define FPROP_ALT_INCLUDED

#include "lib/Measurements/Fermion/fprop.h"
#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
#include "lib/Smear/director_Smear.h"

#include "lib_alt/Solver/asolver.h"

//! Get quark propagator for Fopr with lexical site index: alternative version.

/*!
    This is temporary implementation.
                                        [30 May 2017 H.Matsufuru]
 */

template<typename AFIELD>
class Fprop_alt : public Fprop
{
 public:
  typedef typename AFIELD::real_t real_t;
  //  static const std::string class_name;
  using Fprop::m_vl;
  using Fprop::m_mode;

  virtual ~Fprop_alt() { }

  virtual void set_config(Field *) = 0;

  virtual void invert(Field&, const Field&, int&, double&) = 0;

  virtual void invert_D(Field&, const Field&, int&, double&) = 0;

  virtual void invert_DdagD(Field&, const Field&, int&, double&) = 0;

  virtual void invert(AFIELD&, const AFIELD&, int&, double&) { }

  virtual double flop_count() { return 0.0; }

  virtual void reset_performance() { }

  virtual void get_performance(double& flop_count, double& elapsed_time)
  { }

  virtual void report_performance() { }

  virtual void mult_performance(const std::string mode, const int Nrepeat)
  { }
};
#endif
