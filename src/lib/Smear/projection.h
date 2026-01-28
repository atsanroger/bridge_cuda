/*!
        @file    projection.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef PROJECTION_INCLUDED
#define PROJECTION_INCLUDED

#include "Smear/aprojection.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! base class for projection operator into gauge group.

/*!
  The base class for projection operator into gauge group.

 - original version               [07 Apr 2012 H.Matsufuru]
 - typedef of AProjection<Field>  [07 Feb 2023 H.Matsufuru]
 */

typedef AProjection<Field> Projection;

/*
class Projection
{
 public:
  Projection() {}
  virtual ~Projection() {}

 private:
  // non-copyable
  Projection(const Projection&);
  Projection& operator=(const Projection&);

 public:
  //! projection V = P[alpha, C, U]
  virtual void project(Field_G& v,
                       const double alpha,
                       const Field_G& C, const Field_G& U) = 0;

  //! determination of fields for force calculation
  virtual void force_recursive(Field_G& Xi, Field_G& iTheta,
                               const double alpha, const Field_G& Sigmap,
                               const Field_G& C, const Field_G& U) = 0;

  virtual void set_parameters(const Parameters& param) = 0;

  virtual void get_parameters(Parameters& param) const = 0;


#ifdef USE_FACTORY
 public:
  typedef Projection *(*ProductCreator)();
  typedef Projection *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<Projection, ProductCreator> Factory;
  typedef FactoryTemplate<Projection, ProductCreator_params> Factory_params;

  static Projection *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static Projection *New(const IdentifierType& subtype, const Parameters& params)
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

*/

#endif
