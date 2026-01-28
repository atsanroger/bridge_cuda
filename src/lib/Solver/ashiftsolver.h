/*!
        @file    ashiftsolver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-06-29 22:43:06 #$

        @version $LastChangedRevision: 2527 $
*/

#ifndef ASHIFTSOLVER_INCLUDED
#define ASHIFTSOLVER_INCLUDED

#include "bridge_defs.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
#include "Field/field.h"
#include "Fopr/fopr.h"

//! Shiftsolver class as an abstract base class for multi-shift solvers.
template<typename AFIELD>
class AShiftsolver
{
 public:
  typedef typename AFIELD::real_t real_t;

  AShiftsolver() {}

  virtual ~AShiftsolver() {}

 private:
  // non-copyable
  AShiftsolver(const AShiftsolver&);
  AShiftsolver& operator=(const AShiftsolver&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void solve(std::vector<AFIELD>& solution,
                     const std::vector<real_t>& shift,
                     const AFIELD& source,
                     int& Nconv,
                     real_t& diff) = 0;

  virtual double flop_count() = 0;
};

#endif
