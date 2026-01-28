/*!
        @file    aforceSmear.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef AFORCESMEAR_INCLUDED
#define AFORCESMEAR_INCLUDED

#include "lib/Smear/projection.h"
#include "lib/Field/field_G.h"

#include "lib/IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif


//! Base class for force calculation of smeared operators.

/*!
   Base class for force calculation of smeared operators.
                                     [28 Dec 2011 H.Matsufuru]
   - converted to template base class.
                                     [21 Mar 2023 H.Matsufuru]
 */

template<typename AFIELD>
class AForceSmear
{
 public:

  AForceSmear() {}
  //  : m_vl(CommonParameters::Vlevel()) {}

  virtual ~AForceSmear() {}

 private:
  // non-copyable
  AForceSmear(const AForceSmear&);
  AForceSmear& operator=(const AForceSmear&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  //void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void get_parameters(Parameters&) const = 0;

  virtual void force_udiv(Field_G&, const Field_G&, const Field_G&) {}

  //protected:
  //Bridge::VerboseLevel m_vl;


#ifdef USE_FACTORY
 public:
  typedef AForceSmear *(*ProductCreator)(AProjection<AFIELD> *);
  typedef AForceSmear *(*ProductCreator_params)(AProjection<AFIELD>*,
                                                const Parameters& params);

  typedef FactoryTemplate<AForceSmear, ProductCreator>          Factory;
  typedef FactoryTemplate<AForceSmear, ProductCreator_params>   Factory_params;

  static AForceSmear *New(const IdentifierType& subtype,
                          AProjection<AFIELD> *proj)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)(proj) : 0;
  }

  static AForceSmear *New(const IdentifierType& subtype,
                          AProjection<AFIELD> *proj, const Parameters& params)
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
#endif
