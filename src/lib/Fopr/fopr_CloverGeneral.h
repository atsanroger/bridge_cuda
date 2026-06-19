/*!
        @file    fopr_CloverGeneral.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#ifndef FOPR_CLOVERGENERAL_INCLUDED
#define FOPR_CLOVERGENERAL_INCLUDED

#include "fopr_WilsonGeneral.h"
#include "fopr_CloverTerm_General.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover general fermion operator.

/*!
    This class implements Clover General fermion operator,
    including Clover fermion on anisotropic lattice,
      T. R. Klassen, Nucl. Phys. B509 391 (1998)
      QCD-TARO Collab., Nucl. Phys. Proc. Suppl. 83 411 (2000)
    and relativistic heavy quarks,
      A. X. El-Khadra et al., Phys. Rev. D 55, 3933 (1997)
      S. Aoki et al., Prog. Theor. Phys. 109, 383 (2003)
      N. H. Christ et al., Phys. Rev. D 76, 074505 (2007)
                                     [21 Mar 2015 Y.Namekawa]
 */

class Fopr_CloverGeneral : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_kappa_s;                       //!< spatial hopping parameter
  double m_kappa_t;                       //!< temporal hopping parameter
  double m_nu_s;                          //!< spatial anisotopy parameter
  double m_r_s;                           //!< spatial Wilson parameter
  double m_cSW_s;                         //!< spatial clover coefficient
  double m_cSW_t;                         //!< temporal clover coefficient
  std::vector<int> m_boundary;            //!< boundary conditions
  std::string m_repr;                     //!<  gamma matrix representation
  std::string m_mode;                     //!<  mode of multiplication

  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF; //!< internal parameters

  Fopr_WilsonGeneral *m_fopr_w;           //!< Wilson fermion kernel
  Fopr_CloverTerm_General *m_fopr_csw;    //!< Clover term operator
  const Field_G *m_U;                     //!< gauge configuration (pointer)

  //  Fopr_Clover_imp* m_imp;    //!< pimple prescription
  Field m_v1, m_v2;              //!< working field.

 public:
  DEPRECATED
  Fopr_CloverGeneral()
    : m_vl(CommonParameters::Vlevel())
  {
    init("Dirac");
  }

  DEPRECATED
  Fopr_CloverGeneral(const std::string repr)
    : m_vl(CommonParameters::Vlevel())
  {
    init(repr);
  }

  Fopr_CloverGeneral(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");
    init(repr);
    set_parameters(params);
  }

  ~Fopr_CloverGeneral()
  {
    tidyup();
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa_s, const double kappa_t,
                      const double nu_s, const double r_s,
                      const double cSW_s, const double cSW_t,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
    m_fopr_csw->set_config(U);
  }

  void set_mode(const std::string mode);

  std::string get_mode() const
  {
    return m_mode;
  }

  void mult(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      D(v, f);
    } else if (m_mode == "Ddag") {
      Ddag(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "DDdag") {
      DDdag(v, f);
    } else if (m_mode == "H") {
      H(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void mult_dag(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      Ddag(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "Ddag") {
      D(v, f);
    } else if (m_mode == "H") {
      H(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void DdagD(Field&, const Field&);
  void DDdag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w)
  {
    m_fopr_w->mult_gm5(v, w);
  }

  void mult_up(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_up(mu, v, w);
  }

  void mult_dn(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_dn(mu, v, w);
  }

  void mult_isigma(Field_F&, const Field_F&,
                   const int mu, const int nu);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nex() { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init(const std::string repr);
  void tidyup();

#ifdef USE_FACTORY
 private:
  static Fopr *create_object()
  {
    return new Fopr_CloverGeneral();
  }

  static Fopr *create_object_with_arg(const std::string& repr)
  {
    return new Fopr_CloverGeneral(repr);
  }

  static Fopr *create_object_with_params(const Parameters& params)
  {
    return new Fopr_CloverGeneral(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("CloverGeneral", create_object);
    init &= Fopr::Factory_string::Register("CloverGeneral", create_object_with_arg);
    init &= Fopr::Factory_params::Register("CloverGeneral", create_object_with_params);
    return init;
  }
#endif
};
#endif
