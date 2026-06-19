/*!
        @file    fopr_Clover_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_CLOVER_SF_INCLUDED
#define FOPR_CLOVER_SF_INCLUDED

#include "fopr_Wilson_SF.h"
#include "Measurements/Gauge/staple_SF.h"

#include "Field/shiftField_lex.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover fermion operator.

/*!
    This class implements the clover (improved Wilson) fermion
    operator with SF BC.
    <ul>
    <li>The field strength is calculate when the function set_config() is called.
    <li>Dirac representation only!
    <li>[10 Apr 2012 Y.Taniguchi]
    </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.        [29 Dec 2022 H.Matsufuru]
 */


class Fopr_Clover_SF : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF;
  double m_kappa, m_cSW;
  std::vector<int> m_boundary;
  std::string m_repr;
  std::string m_mode;

  Fopr_Wilson_SF *m_fopr_w;
  const Field_G *m_U;
  ShiftField_lex m_shift;

  Field_G m_Bx, m_By, m_Bz, m_Ex, m_Ey, m_Ez;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

  std::vector<GammaMatrix> m_GM, m_SG;

  std::vector<double> m_phi;    //!< SF boundary condition at t=0
  std::vector<double> m_phipr;  //!< SF boundary condition at t=Nt

  Field m_w1, m_w2;

 public:
  DEPRECATED
  Fopr_Clover_SF()
    : m_vl(CommonParameters::Vlevel())
  {
    init("Dirac");
  }

  Fopr_Clover_SF(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");

    if (repr != "Dirac") {
      vout.crucial(m_vl, "only Dirac representation is supported.\n");
      exit(EXIT_FAILURE);
    }

    init(repr);
    set_parameters(params);
  }

  ~Fopr_Clover_SF() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const double cSW, const std::vector<int> bc,
                      const std::vector<double> phi,
                      const std::vector<double> phipr);

  void get_parameters(Parameters& params) const;

  //! setup configuration (Note that this method is not multi-threaded).
  void set_config(Field *U);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f);

  void mult(Field& v, const Field& f, const std::string mode);

  void mult_dag(Field& v, const Field& f, const std::string mode);

  void DdagD(Field&, const Field&);
  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w)
  { m_fopr_w->mult_gm5(v, w); }

  void mult_isigma(Field_F&, const Field_F&,
                   const int mu, const int nu);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nex() { return 1; }

  //! this returns the number of floating point number operations.
  double flop_count();

 private:
  void init(const std::string repr);
  void tidyup();

  void set_csw();
  void mult_csw(Field&, const Field&);
  void set_fieldstrength(Field_G&, const int, const int);

  void mult_csw_dirac(Field&, const Field&);

  //void mult_csw_chiral(Field&, const Field&);

  void set_boundary_zero(Field&);

  int sg_index(const int mu, const int nu) { return mu * m_Ndim + nu; }


#ifdef USE_FACTORY
 private:

  static Fopr *create_object()
  {
    return new Fopr_Clover_SF();
  }

  static Fopr *create_object_with_params(const Parameters& params)
  {
    return new Fopr_Clover_SF(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Clover_SF", create_object);
    init &= Fopr::Factory_params::Register("Clover_SF", create_object_with_params);
    return init;
  }
#endif
};
#endif
