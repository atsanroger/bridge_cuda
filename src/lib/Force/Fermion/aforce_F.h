/*!
        @file    aforce_F.h
 
        @brief
                                                                               
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
                                                                               
        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$
                                                                               
        @version $LastChangedRevision: 2513 $
*/

#ifndef AFORCE_F_INCLUDED
#define AFORCE_F_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Fopr/afopr.h"
#include "lib/Field/field.h"
#include "lib/Field/field_G.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Tools/director.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "lib/Tools/factory.h"
#endif

/*!
    The base template class of fermion force.
                                   [06 Nov 2018 H.Matsufuru]
*/

template<typename AFIELD>
class AForce_F
{
 protected:
  //Bridge::VerboseLevel m_vl;
  static const std::string class_name;

  Field_G *m_U;   //!< Gauge configuration
  AFIELD m_Ucp;   //!< Gauge configuration (converted if necessary)
  AFIELD m_ut;    //!< working vector for gauge field

 public:
  AForce_F()
    : m_U(0) { init(); }

  virtual ~AForce_F() { tidyup(); }

 private:
  //! non-copyable
  AForce_F(const AForce_F<AFIELD>&);

  //! non-copyable
  AForce_F& operator=(const AForce_F&);

 public:

  //! sets parameters by a Parameter object: to be implemented in a subclass.
  virtual void set_parameters(const Parameters& params){
    vout.crucial("AFopr: set_parameters not implemented.\n");
  }

  virtual void get_parameters(Parameters& params) const {
    vout.crucial("AFopr: get_parameters not implemented.\n");
  }

  //! sets verbose level.
  //void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  //{ m_vl = vl; }

  //! sets the gauge configuration.
  virtual void set_config(Field*) = 0;

  //! in Force, setting the mode is optional when H is nonhermitian.
  virtual void set_mode(const std::string& mode)
  {  /* do nothing if not defined in a subclass. */  }


  //  virtual void force_core(Field&, const Field&);
  virtual void force_core(AFIELD&, const AFIELD&);

  //  virtual void force_core1(Field&, const Field&, const Field&);
  virtual void force_core1(AFIELD&, const AFIELD&, const AFIELD&);

  //  virtual void force_udiv(Field&, const Field&) {}
  virtual void force_udiv(AFIELD&, const AFIELD&) {}

  //  virtual void force_udiv1(Field&, const Field&, const Field&) {}
  virtual void force_udiv1(AFIELD&, const AFIELD&, const AFIELD&) {}

 private:
  //! initializer.
  virtual void init();

  //! finalizer.
  virtual void tidyup();

  virtual void mult_generator(AFIELD&) {}

  virtual void mult_generator(Field_G&) {}


#ifdef USE_FACTORY
 public:
  typedef AForce_F *(*ProductCreator_noarg)();
  typedef AForce_F *(*ProductCreator_string)(const std::string& arg);
  typedef AForce_F *(*ProductCreator_params)(const Parameters& params);
  typedef AForce_F *(*ProductCreator_fopr_params)
                       (AFopr<AFIELD> *fopr, const Parameters& params);
  typedef AForce_F *(*ProductCreator_force_director)
                          (AForce_F<AFIELD> *force_F, Director *director);
  typedef AForce_F *(*ProductCreator_force_forcesmear_director)
                          (AForce_F<AFIELD> *force_F,
                           ForceSmear *forceSmear, Director *director);
  typedef AForce_F *(*ProductCreator_fopr_force)
                       (AFopr<AFIELD> *fopr, AForce_F<AFIELD> *force_F);

  typedef FactoryTemplate<AForce_F, ProductCreator_noarg>  Factory_noarg;
  typedef FactoryTemplate<AForce_F, ProductCreator_string> Factory_string;
  typedef FactoryTemplate<AForce_F, ProductCreator_params> Factory_params;
  typedef FactoryTemplate<AForce_F, ProductCreator_fopr_params>
                                                      Factory_fopr_params;
  typedef FactoryTemplate<AForce_F, ProductCreator_force_director>
                                                   Factory_force_director;
  typedef FactoryTemplate<AForce_F, ProductCreator_force_forcesmear_director>
                                       Factory_force_forcesmear_director;
  typedef FactoryTemplate<AForce_F, ProductCreator_fopr_force>
                                                       Factory_fopr_force;

  static AForce_F *New(const IdentifierType& subtype)
  {
    ProductCreator_noarg p = Factory_noarg::Find(subtype);
    return p ? (*p)() : 0;
  }

  static AForce_F *New(const IdentifierType& subtype, const std::string& arg)
  {
    ProductCreator_string p = Factory_string::Find(subtype);
    return p ? (*p)(arg) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype,
                       AFopr<AFIELD> *fopr, const Parameters& params)
  {
    ProductCreator_fopr_params p = Factory_fopr_params::Find(subtype);
    return p ? (*p)(fopr, params) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype,
                       AForce_F<AFIELD> *force_F, Director *director)
  {
    ProductCreator_force_director p = Factory_force_director::Find(subtype);
    return p ? (*p)(force_F, director) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype,
                       AForce_F<AFIELD> *force_F, ForceSmear* forcesmear,
		       Director *director)
  {
    ProductCreator_force_forcesmear_director p
                    = Factory_force_forcesmear_director::Find(subtype);
    return p ? (*p)(force_F, forcesmear, director) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype,
                       AFopr<AFIELD> *fopr, AForce_F<AFIELD> *force_F)
  {
    ProductCreator_fopr_force p = Factory_fopr_force::Find(subtype);
    return p ? (*p)(fopr, force_F) : 0;
  }

  /*
  static AForce_F *New(const IdentifierType& subtype, AFopr *fopr)
  {
    ProductCreator_fopr p = Factory_fopr::Find(subtype);
    return p ? (*p)(fopr) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype, unique_ptr<AFopr>& fopr)
  {
    ProductCreator_fopr p = Factory_fopr::Find(subtype);
    return p ? (*p)(fopr.get()) : 0;
  }

  static AForce_F *New(const IdentifierType& subtype, unique_ptr<AFopr>& fopr,
		    unique_ptr<Director>& director)
  {
    ProductCreator_fopr_director p = Factory_fopr_director::Find(subtype);
    return p ? (*p)(fopr.get(), director.get()) : 0;
  }
  */


#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif  // USE_FACTORY_AUTOREGISTER

#endif  // USE_FACTORY

};
#endif  // AFORCE_F_INCLUDED

