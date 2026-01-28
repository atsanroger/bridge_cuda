/*!
        @file    afopr_Smeared.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#ifndef AFOPR_SMEARED_INCLUDED
#define AFOPR_SMEARED_INCLUDED

#include "lib/Fopr/afopr.h"
#include "lib/Smear/director_Smear.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! smeared fermion operator: alternative version.

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
template<typename AFIELD>
class AFopr_Smeared : public AFopr<AFIELD>
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  AFopr<AFIELD> *m_fopr;
  Director_Smear *m_dr_smear;

 public:
  //! constructor requires Fopr and Director_Smear objects
  AFopr_Smeared(AFopr<AFIELD> *fopr, Director_Smear *dr_smear)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_dr_smear(dr_smear) {}

  AFopr_Smeared(AFopr<AFIELD> *fopr, Director_Smear *dr_smear, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_fopr(fopr), m_dr_smear(dr_smear)
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  //! set pointer to original thin link variable
  void set_config(Field *U);

  //! multiply smeared fermion operator
  void mult(AFIELD& v, const AFIELD& f)
  { m_fopr->mult(v, f); }

  //! multiply smeared fermion operator
  void mult_dag(AFIELD& v, const AFIELD& f)
  { m_fopr->mult_dag(v, f); }

  //! multiply smeared fermion operator
  void mult(AFIELD& v, const AFIELD& f, std::string mode)
  { m_fopr->mult(v, f, mode); }

  //! multiply smeared fermion operator
  void mult_dag(AFIELD& v, const AFIELD& f, std::string mode)
  { m_fopr->mult_dag(v, f, mode); }

  //! multiply gamma_5 matrix.
  void mult_gm5(AFIELD& v, const AFIELD& f)
  { m_fopr->mult_gm5(v, f); }

  //! set the mode of fermion operator.
  void set_mode(const std::string mode)
  { m_fopr->set_mode(mode); }

  //! requirement of spinor field conversion.
  bool needs_convert()
  { return m_fopr->needs_convert(); }

  //! convert of spinor field.
  void convert(AFIELD& v, const Field& w)
  { m_fopr->convert(v, w); }

  //! reverse of spinor field.
  void reverse(Field& v, const AFIELD& w)
  { m_fopr->reverse(v, w); }

  std::string get_mode() const
  { return m_fopr->get_mode(); }

  void mult_up(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr->mult_up(mu, v, w); }

  void mult_dn(const int mu, AFIELD& v, const AFIELD& w)
  { m_fopr->mult_dn(mu, v, w); }

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin()  { return m_fopr->field_nin(); }
  int field_nex()  { return m_fopr->field_nex(); }

  //! returns floating operation counts.
  double flop_count();

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object(AFopr<AFIELD> *fopr, Director *director)
  {
    return new AFopr_Smeared(fopr, dynamic_cast<Director_Smear *>(director));
  }

  static AFopr<AFIELD> *create_object_with_params(AFopr<AFIELD> *fopr, Director *director, const Parameters& params)
  {
    return new AFopr_Smeared(fopr, dynamic_cast<Director_Smear *>(director), params);
  }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_fopr_director::Register("Smeared", create_object);
    // init &= AFopr<AFIELD>::Factory_fopr_director_params::Register("Smeared", create_object_with_params);
    return init;
  }
#endif
};
#endif
