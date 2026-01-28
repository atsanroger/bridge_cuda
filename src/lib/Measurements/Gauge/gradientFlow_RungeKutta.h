/*!
        @file    gradientFlow_RungeKutta.h

        @brief

        @author  Yusuke Namekawa (namekawa)

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef GRADIENTFLOW_RUNGEKUTTA_INCLUDED
#define GRADIENTFLOW_RUNGEKUTTA_INCLUDED

#include "Action/action.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! GradientFlow_RungeKutta construction.

/*!
    This class implements Runge-Kutta method
    for GradientFlow.
                              [10 Oct 2014 Y.Namekawa]
 */

class GradientFlow_RungeKutta
{
 public:
  GradientFlow_RungeKutta(Action *action, const int Nprec, const Bridge::VerboseLevel vl)
  {}

  virtual ~GradientFlow_RungeKutta() {}

  virtual void flow(double& t, double& Estep, Field_G& U) = 0;
  virtual int Norder_RK() const = 0;
};
#endif
