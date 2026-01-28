/*!
        @file    fopr_Smeared.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FOPR_SMEARED_INCLUDED
#define FOPR_SMEARED_INCLUDED

#include "fopr.h"
#include "Smear/director_Smear.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! smeared fermion operator.

/*!
    This class construct a smeared fermon operator for a given
    base fermion operator together with smearing director.
    Both of them must be constructed beforehand outside this
    class and given to the constructor.
    Smearing of link configuration is triggered by call of
    set_config(), which calls set_config() of the smearing
    director and then gets the pointer to the smeared
    config. to set it as the config. of base fermion operator.
    When mult() is called, this class just call the mult()
    of base class.
                                  [24 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                  [21 Mar 2015 Y.Namekawa]
 */

class Fopr_Smeared : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Fopr *m_fopr;
  Director_Smear *m_dr_smear;

 public:
  //! constructor requires Fopr and Director_Smear objects
  Fopr_Smeared(Fopr *fopr, Director_Smear *dr_smear)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_dr_smear(dr_smear) {}

  Fopr_Smeared(Fopr *fopr, Director_Smear *dr_smear, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_dr_smear(dr_smear)
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  //! set pointer to original thin link variable
  void set_config(Field *U);

  //! multiply smeared fermion operator
  void mult(Field& v, const Field& f)
  {
    m_fopr->mult(v, f);
  }

  //! multiply smeared fermion operator
  void mult_dag(Field& v, const Field& f)
  {
    m_fopr->mult_dag(v, f);
  }

  //! set the mode of fermion operator
  void set_mode(const std::string mode)
  {
    m_fopr->set_mode(mode);
  }

  std::string get_mode() const
  {
    return m_fopr->get_mode();
  }

  void mult_up(const int mu, Field& v, const Field& w)
  {
    m_fopr->mult_up(mu, v, w);
  }

  void mult_dn(const int mu, Field& v, const Field& w)
  {
    m_fopr->mult_dn(mu, v, w);
  }

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin() { return m_fopr->field_nin(); }
  int field_nex() { return m_fopr->field_nex(); }

  //! this returns the number of floating point operations.
  double flop_count();

#ifdef USE_FACTORY
 private:
  static Fopr *create_object(Fopr *fopr, Director *director)
  {
    return new Fopr_Smeared(fopr, dynamic_cast<Director_Smear *>(director));
  }

  static Fopr *create_object_with_params(Fopr *fopr, Director *director, const Parameters& params)
  {
    return new Fopr_Smeared(fopr, dynamic_cast<Director_Smear *>(director), params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_fopr_director::Register("Smeared", create_object);
    // init &= Fopr::Factory_fopr_director_params::Register("Smeared", create_object_with_params);
    return init;
  }
#endif
};
#endif
