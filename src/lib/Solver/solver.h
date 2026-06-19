/*!
        @file    solver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

#include "Fopr/fopr.h"

#include "ResourceManager/threadManager.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class for linear solver class family.

/*!
                               [28 Dec 2011 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Add restart.               [22 Feb 2016 Y.Namekawa]
    Add flop_count.            [ 8 Aug 2016 Y.Namekawa]
    Add use_init_guess.       [ 7 Jul 2017 Y.Namekawa]
 */

class Solver
{
 public:
  Solver() {}

  virtual ~Solver() {}

 private:
  Solver(const Solver&);
  Solver& operator=(const Solver&);

 public:
  virtual void set_parameters(const Parameters& params) = 0;

  virtual void set_parameters(const int Niter, const int Nrestart, const double Stop_cond) = 0;

  // virtual void set_parameters(const int Niter, const int Nrestart, const double Stop_cond, const bool use_init_guess) = 0;

  virtual void get_parameters(Parameters& params) const = 0;

  virtual void solve(Field& solution, const Field& source, int& Nconv, double& diff) = 0;

  virtual Fopr *get_fopr() = 0;

  virtual double flop_count() = 0;


#ifdef USE_FACTORY
 public:
  typedef Solver *(*ProductCreator)(Fopr *);
  typedef Solver *(*ProductCreator_params)(Fopr *, const Parameters&);

  typedef FactoryTemplate<Solver, ProductCreator>          Factory;
  typedef FactoryTemplate<Solver, ProductCreator_params>   Factory_params;

  static Solver *New(const IdentifierType& subtype, Fopr *fopr)
  {
    ProductCreator p = Factory::Find(subtype);

    return p ? (*p)(fopr) : 0;
  }

  static Solver *New(const IdentifierType& subtype, Fopr *fopr, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);

    return p ? (*p)(fopr, params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
