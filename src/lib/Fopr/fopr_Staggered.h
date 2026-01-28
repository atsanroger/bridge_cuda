/*!
        @file    fopr_Staggered.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_STAGGERED_INCLUDED
#define FOPR_STAGGERED_INCLUDED

#include "Fopr/fopr.h"
//#include "Field/field_F_1spinor.h"

//#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Staggered fermion operator.

/*!
    This class is a complexified standard staggered fermion operator.
    The implementation is based on the even-odd staggered version.
                                        [19 Sep 2015 H.Matsufuru]
 */

class Fopr_Staggered : public Fopr
{
 public:
  static const std::string class_name;

 private:
  // input parameters
  double m_mq;                 //!< quark mass.
  std::vector<int> m_boundary; //!< boundary conditions.
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::string m_mode;          //!< mult mode

  // internal data members
  int m_Nc, m_Nvc, m_Ndf;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol, m_Ndim;

  Field m_stg_phase; //!< staggered phase
  Field m_parity;    //!< site parity field

  Field_G m_U;       //!< gauge field multiplied by staggered phase.

  Field m_v1;        //!< working field
  Field m_v2;        //!< working field

  //! communication buffer
  double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
  double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
  double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
  double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

 public:
  //! constructor.
  Fopr_Staggered(const Parameters& params) { init(params); }

  DEPRECATED
  Fopr_Staggered() : Fopr() { init(); }

  //! destructor.
  ~Fopr_Staggered() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(Field&, const Field&);

  void mult_dag(Field&, const Field&);

  void mult_gm5(Field&, const Field&);

  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void DdagD(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field&);

  void mult_staggered_phase(Field&, int mu);

  void normalize_fprop(Field& v) { scal(v, 1.0 / m_mq); }

  void normalize_fopr_(Field& v) { scal(v, m_mq); }

  int field_nin() { return 2 * m_Nc; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return 1; }

  //! returns the number of floating point operations.
  double flop_count() { return flop_count(m_mode); }

  //! returns the flop count for specified mode.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  void init();

  void init(const Parameters& params);

  void setup();

  //! final clean-up.
  void tidyup();

  void set_staggered_phase();

  void set_config_omp(Field *U);

  void set_config_impl(Field *U);

  void mult_xp(Field&, const Field&);
  void mult_xm(Field&, const Field&);
  void mult_yp(Field&, const Field&);
  void mult_ym(Field&, const Field&);
  void mult_zp(Field&, const Field&);
  void mult_zm(Field&, const Field&);
  void mult_tp(Field&, const Field&);
  void mult_tm(Field&, const Field&);

#ifdef USE_FACTORY
 private:
  static Fopr *create_object() { return new Fopr_Staggered(); }

  static Fopr *create_object_with_params(const Parameters& params)
  { return new Fopr_Staggered(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Staggered", create_object);
    init &= Fopr::Factory_params::Register("Staggered",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
