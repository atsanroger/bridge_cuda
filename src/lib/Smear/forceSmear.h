/*!
        @file    forceSmear.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FORCESMEAR_INCLUDED
#define FORCESMEAR_INCLUDED

#include "lib/Smear/aforceSmear.h"
#include "lib/Field/field.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class for force calculation of smeared operators.

/*!
                                     [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

typedef AForceSmear<Field> ForceSmear;
/*
class ForceSmear
{
 public:

  ForceSmear()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~ForceSmear() {}

 private:
  // non-copyable
  ForceSmear(const ForceSmear&);
  ForceSmear& operator=(const ForceSmear&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void get_parameters(Parameters&) const = 0;

  virtual void force_udiv(Field_G&, const Field_G&, const Field_G&) {}

 protected:
  Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
 public:
  typedef ForceSmear *(*ProductCreator)(Projection *);
  typedef ForceSmear *(*ProductCreator_params)(Projection *, const Parameters& params);

  typedef FactoryTemplate<ForceSmear, ProductCreator>          Factory;
  typedef FactoryTemplate<ForceSmear, ProductCreator_params>   Factory_params;

  static ForceSmear *New(const IdentifierType& subtype, Projection *proj)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)(proj) : 0;
  }

  static ForceSmear *New(const IdentifierType& subtype, Projection *proj, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(proj, params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};

*/

#endif
