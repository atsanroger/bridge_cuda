/*!
        @file    aeigensolver.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef AEIGENSOLVER_INCLUDED
#define AEIGENSOLVER_INCLUDED

#include "bridge_defs.h"
#include "Parameters/commonParameters.h"
#include "Parameters/parameters.h"
#include "complexTraits.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "Tools/factory.h"
#include "Tools/director.h"
#endif

//! Eigensolver class for abstract base class of eigen solvers.

/**
   Eigensolver class provides an abstract base class for solvers
   of eigenvalues and eigenvectors.

   This class was changed to template class so as to be available
   with general gield and fermion operator classes.
   The original claas was Eigensolver.
                                       [23 Apr 2018 H.Matsufuru]
 */

template<typename FIELD, typename FOPR>
class AEigensolver
{
 public:
  typedef typename FIELD::real_t                      real_t;
  typedef typename ComplexTraits<real_t>::complex_t   complex_t;

  AEigensolver() {}

  virtual ~AEigensolver() {}

 private:
  // non-copyable
  AEigensolver(const AEigensolver<FIELD, FOPR>&);
  AEigensolver& operator=(const AEigensolver<FIELD, FOPR>&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void solve(std::vector<real_t>& TDa,
                     std::vector<FIELD>& vk,
                     int& Nsbt, int& Nconv,
                     const FIELD& b)
  {
    vout.crucial("AEigensolver: real solve not implemented.\n");
    exit(EXIT_FAILURE);
  }

  //! complex version of solve.
  virtual void solve(std::vector<complex_t>& TDa,
                     std::vector<FIELD>& vk,
                     int& Nsbt, int& Nconv,
                     const FIELD& b)
  {
    vout.crucial("AEigensolver: complex solve not implemented.\n");
    exit(EXIT_FAILURE);
  }

#ifdef USE_FACTORY
 public:
  typedef AEigensolver *(*ProductCreator)(FOPR *);
  typedef AEigensolver *(*ProductCreator_params)(FOPR *, const Parameters& params);

  typedef FactoryTemplate<AEigensolver, ProductCreator>          Factory_fopr;
  typedef FactoryTemplate<AEigensolver, ProductCreator_params>   Factory_fopr_params;

  static AEigensolver *New(const IdentifierType& subtype, FOPR *fopr)
  {
    ProductCreator p = Factory_fopr::Find(subtype);
    return p ? (*p)(fopr) : 0;
  }

  static AEigensolver *New(const IdentifierType& subtype, FOPR *fopr, const Parameters& params)
  {
    ProductCreator_params p = Factory_fopr_params::Find(subtype);
    return p ? (*p)(fopr, params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif  // USE_FACTORY
};

#endif
