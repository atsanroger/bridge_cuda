/*!
        @file    aprojection.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef APROJECTION_INCLUDED
#define APROJECTION_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif

//! Base template class for projection operator into gauge group.

/*!
  This template class provides the base class of projection into
  SU(N) gauge field, which is assumed to be combined with smearing
  operators.
 
  - Original version                    [07 Apr 2012 H.Matsufuru]
  - template base class                 [07 Feb 2023 H.Matsufuru]
 */
template<typename AFIELD>
class AProjection
{
 public:
  AProjection() {}
  virtual ~AProjection() {}

 private:
  // non-copyable
  AProjection(const AProjection&);
  AProjection& operator=(const AProjection&);

 public:
  //! projection V = P[alpha, C, U]
  virtual void project(Field_G& v,
                       const double alpha,
                       const Field_G& C, const Field_G& U) = 0;

  virtual void project(AFIELD& v,
                       const double alpha,
                       const AFIELD& C, const AFIELD& U)
  {  };

  //! determination of fields for force calculation
  virtual void force_recursive(Field_G& Xi, Field_G& iTheta,
                               const double alpha, const Field_G& Sigmap,
                               const Field_G& C, const Field_G& U) = 0;

  virtual void set_parameters(const Parameters& param) = 0;

  virtual void get_parameters(Parameters& param) const = 0;


#ifdef USE_FACTORY
 public:
  typedef AProjection *(*ProductCreator)();
  typedef AProjection *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<AProjection, ProductCreator> Factory;
  typedef FactoryTemplate<AProjection, ProductCreator_params> Factory_params;

  static AProjection *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static AProjection *New(const IdentifierType& subtype,
                          const Parameters& params)
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
