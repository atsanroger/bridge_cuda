/*!
        @file    force_F_Staggered_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#ifndef FORCE_F_STAGGERED_EO_INCLUDED
#define FORCE_F_STAGGERED_EO_INCLUDED

#include "force_F.h"
#include "Fopr/fopr_Staggered_eo.h"
#include "Field/field_F_1spinor.h"

//! Force calculation for even-odd staggered fermion operator.

/*!
    This is a temporary implementation.
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
 */


class Force_F_Staggered_eo : public Force
{
 public:
  static const std::string class_name;

 private:
  // input parameters
  double m_mq;
  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  // internal data members
  int m_Nvol, m_Ndim;

  Fopr_Staggered_eo *m_fopr_ks;
  Field_F_1spinor m_psf;
  Index_eo m_index_eo;
  ShiftField_eo m_shift_eo;

  Field_G *m_Ueo;

 public:
  Force_F_Staggered_eo() { init(); }

  Force_F_Staggered_eo(const Parameters& params)
  {
    init();
    set_parameters(params);
  }

  ~Force_F_Staggered_eo() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const double mq, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_core(Field& force, const Field& eta);
  void force_udiv(Field& force, const Field& eta);

  void force_core1(Field&, const Field&, const Field&);  // dummy entry
  void force_udiv1(Field&, const Field&, const Field&);  // dummy entry

  void force_udiv1(Field_G& force, const Field_F_1spinor& zeta, const Field_F_1spinor& eta, int ieo);

 private:
  //! initial setup.
  void init();

  //! final clean-up.
  void tidyup();
};
#endif
