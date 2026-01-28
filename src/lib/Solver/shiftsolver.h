/*!
        @file    shiftsolver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SHIFTSOLVER_INCLUDED
#define SHIFTSOLVER_INCLUDED

/*
#include "defs.h"
#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
#include "Field/field.h"
#include "Fopr/fopr.h"
*/

//! Shiftsolver class as an abstract base class for multi-shift solvers.

#include "Solver/ashiftsolver.h"
#include "Field/field.h"

typedef AShiftsolver<Field> Shiftsolver;

/*
class Shiftsolver
{
 public:
  Shiftsolver() {}

  virtual ~Shiftsolver() {}

 private:
  // non-copyable
  Shiftsolver(const Shiftsolver&);
  Shiftsolver& operator=(const Shiftsolver&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void solve(
    std::vector<Field>& solution,
    const std::vector<double>& shift,
    const Field& source,
    int& Nconv, double& diff) = 0;

};
*/

#endif
