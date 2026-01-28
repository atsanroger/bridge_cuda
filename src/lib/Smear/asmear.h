/*!
        @file    asmear.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef ASMEAR_INCLUDED
#define ASMEAR_INCLUDED

#include "lib/Smear/aprojection.h"
#include "lib/IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif


//! Base template class for smearing of link variables.

/*!
   This template class provide the base class of smearing.

   - original version              [28 Dec 2011 H.Matsufuru]
   - template base class           [20 Feb 2023 H.Matsufuru]
 */

template<typename AFIELD>
class ASmear
{
 public:
  ASmear() {}
  virtual ~ASmear() {}

 private:
  ASmear(const ASmear&);
  ASmear& operator=(const ASmear&);

 public:
  virtual void smear(Field_G&, const Field_G&) = 0;

  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

#ifdef USE_FACTORY
 public:
  typedef ASmear *(*ProductCreator)(AProjection<AFIELD> *);
  typedef ASmear *(*ProductCreator_params)(AProjection<AFIELD> *, const Parameters&);

  typedef FactoryTemplate<ASmear, ProductCreator> Factory;
  typedef FactoryTemplate<ASmear, ProductCreator_params> Factory_params;

  static ASmear *New(const IdentifierType& subtype, AProjection<AFIELD> *proj)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)(proj) : 0;
  }

  static ASmear *New(const IdentifierType& subtype,
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
