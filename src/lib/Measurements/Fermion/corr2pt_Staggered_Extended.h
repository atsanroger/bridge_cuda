/*!
        @file    corr2pt_Staggered_Extended.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef CORR2PT_STAGGERED_EXTENDED_INCLUDED
#define CORR2PT_STAGGERED_EXTENDED_INCLUDED

#include "bridge_defs.h"
#include "Parameters/commonParameters.h"
#include "Field/index_lex.h"
#include "Field/field_F_1spinor.h"
#include "IO/bridgeIO.h"

//! Two-point correlator for staggered fermions.

/*!
   This class implements the staggered hadron correlators
   with non-local sink operators.
   For mesons, quark and antiquark are contracted at
   different points within spatial 2x2x2 cube.
   No phase factor is multiplied for present implementation.
   As for baryons, only the nucleon correlator is implemented
   so far, while Delta-baryon is possible to constructed.
   Cf. Gupta et al., Phys. Rev. D 43, 2003 (1991).
                                    [16 Oct 2020 H.Matsufuru]
   - PS meson correlator is written in the first version.
                                           [28 Dec 2011 H.M.]
   - Baryon correlator is added.           [13 Aug 2018 H.M.]
 */

class Corr2pt_Staggered_Extended
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Index_lex m_index;
  std::vector<int> m_epsilon_index;
  //!< index of totally antisymmetric tensor

  std::string m_filename_output;

 public:
  Corr2pt_Staggered_Extended()
    : m_vl(CommonParameters::Vlevel())
  { init(); }

 private:
  // non-copyable
  Corr2pt_Staggered_Extended(const Corr2pt_Staggered_Extended&);
  Corr2pt_Staggered_Extended& operator=(const Corr2pt_Staggered_Extended&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl)
  { m_vl = vl; }

  void set_output_file(const std::string filename)
  { m_filename_output = filename; }

  void meson_all(std::vector<double>& meson,
                 const std::vector<Field_F_1spinor>& sq1,
                 const std::vector<Field_F_1spinor>& sq2);

  void baryon_all(std::vector<double>& baryon,
                  const std::vector<Field_F_1spinor>& sq1,
                  const std::vector<Field_F_1spinor>& sq2);

  void meson_correlator(std::vector<dcomplex>& corr,
                        const std::vector<int>& shft,
                        const std::vector<Field_F_1spinor>& sq1,
                        const std::vector<Field_F_1spinor>& sq2,
                        const int isrc1, const int isrc2);

  void nucleon_correlator(std::vector<dcomplex>& corr,
                          const std::vector<Field_F_1spinor>& sq1,
                          const std::vector<Field_F_1spinor>& sq2,
                          const int isrc1, const int isrc2);

 private:

  //! initial setup.
  void init();

  //! setup of totally antisymmetric tensor for Nc = 3 case.
  void set_asymtensor_Nc3();

  //! index transformation of totally antisymmetric tensor: index.
  int epsilon_index(int i, int n)
  { return m_epsilon_index[i + 3 * n]; }

  double epsilon_value(int n)
  { return 1.0 - 2.0 * (n / 3); }

  //! color-singlet contraction of three fields in 3D space.
  void contract_at_t(dcomplex& corr, const Field_F_1spinor& v1,
                     const Field_F_1spinor& v2,
                     const Field_F_1spinor& v3,
                     const int time);
};
#endif
