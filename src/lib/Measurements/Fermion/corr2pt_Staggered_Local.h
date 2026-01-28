/*!
        @file    corr2pt_Staggered_Local.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef CORR2PT_STAGGERED_INCLUDED
#define CORR2PT_STAGGERED_INCLUDED

#include "bridge_defs.h"
#include "Parameters/commonParameters.h"
#include "Field/index_lex.h"
#include "Field/field_F_1spinor.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

//! Two-point correlator for staggered fermions.

/*!
   This class implements the local staggered hadron correlators.
   The local correlators mean that the sink operator is composed
   of spatially local contraction of quark propagators.
   For meson, this implies that the PS, SC, VT, AV correlators
   (this notation follows Appendix of Bowler et al., Nucl. Phys.
   B296 (1988) 732) can be described. For baryon, only nucleon
   correlator can be locally composed.
                                        [15 Oct 2020 H.Matsufuru]
   Modified so that odd values of Nx, Ny, Nz are allowed.
                                        [14 Jan 2023 H.Matsufuru]
 */

class Corr2pt_Staggered_Local
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Index_lex m_index;
  std::vector<int> m_epsilon_index;

  std::string m_filename_output;

 public:
  Corr2pt_Staggered_Local()
    : m_vl(CommonParameters::Vlevel()) { init(); }

 private:

  // non-copyable
  Corr2pt_Staggered_Local(const Corr2pt_Staggered_Local&);
  Corr2pt_Staggered_Local& operator=(const Corr2pt_Staggered_Local&);

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
                        const std::vector<double>& prty,
                        const std::vector<Field_F_1spinor>& sq1,
                        const std::vector<Field_F_1spinor>& sq2);

  void nucleon_correlator(std::vector<dcomplex>& corr,
                          const std::vector<Field_F_1spinor>& sq1,
                          const std::vector<Field_F_1spinor>& sq2);

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
