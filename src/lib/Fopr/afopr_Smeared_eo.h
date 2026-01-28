/*!
        @file    afopr_Smeared_eo.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/

#ifndef AFOPR_SMEARED_EO_INCLUDED
#define AFOPR_SMEARED_EO_INCLUDED

#include "Fopr/afopr_eo.h"
#include "Smear/director_Smear.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! smeared fermion operator with even-odd preconditioning.

/*!
    This class constructs a smeared fermon operator for a given
    base fermion operator together with smearing director.
    Both of them must be constructed beforehand outside this
    class and given to the constructor.
    Smearing of link configuration is triggered by call of
    set_config(), which calls set_config() of the smearing
    director and then gets the pointer to the smeared
    config. to set it as the config. of base fermion operator.
    When mult() is called, this class just call the mult()
    of base class.
                                  [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                  [21 Mar 2015 Y.Namekawa]
 */

template<typename AFIELD>
class AFopr_Smeared_eo : public AFopr_eo<AFIELD>
{
 public:
  static const std::string class_name;
  using AFopr_eo<AFIELD>::m_vl;

 private:
  Bridge::VerboseLevel m_vl;

  AFopr_eo<AFIELD> *m_fopr_eo;
  Director_Smear *m_dr_smear;

 public:
  //! constructor requires Fopr and Director_Smear objects
  AFopr_Smeared_eo(AFopr_eo<AFIELD> *fopr_eo, Director_Smear *dr_smear)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr_eo  = fopr_eo;
    m_dr_smear = dr_smear;
  }

  AFopr_Smeared_eo(AFopr_eo<AFIELD> *fopr_eo, Director_Smear *dr_smear, const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    m_fopr_eo  = fopr_eo;
    m_dr_smear = dr_smear;

    set_parameters(params);
  }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  void preProp(AFIELD& Be, AFIELD& bo, const AFIELD& b)
  { m_fopr_eo->preProp(Be, bo, b); }

  void postProp(AFIELD& x, const AFIELD& xe, const AFIELD& bo)
  { m_fopr_eo->postProp(x, xe, bo); }

  //! set pointer to original thin link variable
  void set_config(Field *U);

  //! multiply smeared fermion operator
  void mult(AFIELD& v, const AFIELD& f)
  { m_fopr_eo->mult(v, f); }

  //! multiply smeared fermion operator
  void mult_dag(AFIELD& v, const AFIELD& f)
  { m_fopr_eo->mult_dag(v, f); }

  //! set the mode of fermion operator
  void set_mode(const std::string mode)
  { m_fopr_eo->set_mode(mode); }

  std::string get_mode() const
  { return m_fopr_eo->get_mode(); }

  void mult_up(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr_eo->mult_up(mu, v, w); }

  void mult_dn(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr_eo->mult_dn(mu, v, w); }

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr->convert(v, w); }

  //! reverses a Field object into other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr->reverse(v, w); }

  int field_nvol() { return m_fopr_eo->field_nvol(); }
  int field_nin() { return m_fopr_eo->field_nin(); }
  int field_nex() { return m_fopr_eo->field_nex(); }

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object(AFopr<AFIELD> *fopr, Director *director)
  {
    return new AFopr_Smeared_eo<AFIELD>(
      dynamic_cast<AFopr_eo<AFIELD> *>(fopr),
      dynamic_cast<Director_Smear *>(director));
  }

  static AFopr<AFIELD> *create_object_with_params(AFopr<AFIELD> *fopr, Director *director, const Parameters& params)
  {
    return new AFopr_Smeared_eo<AFIELD>(
      dynamic_cast<AFopr_eo<AFIELD> *>(fopr),
      dynamic_cast<Director_Smear *>(director), params);
  }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_fopr_director::Register("Smeared_eo", create_object);
    // init &= AFopr<AFIELD>::Factory_fopr_director_params::Register("Smeared_eo", create_object_with_params);
    return init;
  }
#endif
};
#endif
