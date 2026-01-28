/*!
        @file    action.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef ACTION_INCLUDED
#define ACTION_INCLUDED

#include "Field/field_G.h"
#include "Tools/randomNumbers.h"
#include "Parameters/parameters.h"
#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class of HMC action class family.

/*!
   This class defines interface of Action-type classes.
                                   [28 Dec 2011 H.Matsufuru]
   Factory is introduced.          [21 Mar 2015 Y.Namekawa]
 */

class Action
{
 public:

  Action() {}

  virtual ~Action() {}

 private:
  // non-copyable
  Action(const Action&);
  Action& operator=(const Action&);

 public:
  virtual void set_parameters(const Parameters& param) = 0;

  virtual void get_parameters(Parameters& param) const = 0;

  //! setting pointer to the gauge configuration.
  virtual void set_config(Field *U) = 0;

  //! Langevis step.
  virtual double langevin(RandomNumbers *) = 0;

  //! calculate Hamiltonian of this action term.
  virtual double calcH() = 0;

  //! returns force for molcular dynamical update of conjugate momenta.
  //virtual const Field force() = 0;
  virtual void force(Field&) = 0;

  virtual void force(Field& v, Field& U)
  {
    set_config(&U);
    force(v);
  }

#ifdef USE_FACTORY
 public:
  typedef Action *(*ProductCreator)();
  typedef Action *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<Action, ProductCreator>          Factory;
  typedef FactoryTemplate<Action, ProductCreator_params>   Factory_params;

  static Action *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static Action *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
