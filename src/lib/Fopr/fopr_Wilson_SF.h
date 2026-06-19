/*!
        @file    fopr_Wilson_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef FOPR_WILSON_SF_INCLUDED
#define FOPR_WILSON_SF_INCLUDED

#include "fopr_Wilson.h"

//#include "Field/field_F_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator with SF BC.

/*!
    This class implements the Wilson fermion operator with the SF BC.
    <ul>
    <li>The SF BC for the fermion field is:
    <ul>
    <li>The fermion field at t=0 is always set to zero (Dirichlet BC).
    <li>The field at t=1,...,Lt-1 is active.
    </ul>
    <li>Implemented by delegation of the Fopr_Wilson class like Fopr_Clover.
    <li>The modification is only in Fopr_Wilson_SF::D to set the boundary fermion field to zero before and after multiplication of the standard Wilson Dirac operator.
    <ul>
    <li>By this manipulation the SF BC is introduced in the fermion field.
    </ul>
    <li>A private function set_boundary_zero(Field&) is introduced for this bounadry manipulation.
    <li>A few private members are added for this function.
    <li> [04 Apr 2012 Y.Taniguchi]
    </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
    Multi-threading appled.         [29 Dec 2022 H.Matsufuru]
 */

class Fopr_Wilson_SF : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF;
  double m_kappa;
  std::vector<int> m_boundary;
  std::string m_repr;

  std::string m_mode;

  Fopr_Wilson *m_fopr_w;
  const Field_G *m_U;

  Field m_w1, m_w2;

 public:

  Fopr_Wilson_SF() { init(); }

  Fopr_Wilson_SF(const Parameters& params)
  {
    init();
    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  ~Fopr_Wilson_SF() { tidyup(); }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f);

  void mult(Field& v, const Field& f, const std::string mode);

  void mult_dag(Field& v, const Field& f, const std::string mode);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void DdagD(Field&, const Field&);
  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w);

  void mult_gm5p(const int mu, Field& v, const Field& w);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nex() { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init();

  void tidyup(); // { delete m_fopr_w; }

  //! A function to set the fermion field to zero at the t=0 boundary.
  void set_boundary_zero(Field&);

#ifdef USE_FACTORY
 private:
  static Fopr *create_object()
  {
    return new Fopr_Wilson_SF();
  }

  static Fopr *create_object_with_params(const Parameters& params)
  {
    return new Fopr_Wilson_SF(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Wilson_SF", create_object);
    init &= Fopr::Factory_params::Register("Wilson_SF",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
