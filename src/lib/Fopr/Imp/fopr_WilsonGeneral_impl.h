/*!
        @file    fopr_WilsonGeneral_impl.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_GENERAL_IMPL_IMP_INCLUDED
#define FOPR_WILSON_GENERAL_IMPL_IMP_INCLUDED

#include "Fopr/fopr.h"

#include "Field/shiftField_lex.h"
#include "Tools/gammaMatrixSet.h"

#include "ResourceManager/threadManager.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson_General fermion operator.

/*!
    This class implements the Wilson_General fermion operator,
    including Wilson fermion on anisotropic lattice, and
    relativistic heavy quarks.
                                    [21 Mar 2015 Y.Namekawa]
    Multi-threading is modified.   [30 Dec 2022 H.Matsufuru]
 */

namespace Imp {
  class Fopr_WilsonGeneral : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    // input parameters
    double m_kappa_s, m_kappa_t;
    double m_nu_s, m_r_s;
    std::vector<int> m_boundary; //!< boundary condition
    std::string m_repr;          //!< gamma-matrix representation
    Bridge::VerboseLevel m_vl;   //!< verbose level

    std::string m_mode;          //!< mult mode

    // internal data members
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nx, m_Ny, m_Nz, m_Nt;
    int m_Nvol, m_Ndim;

    const Field_G *m_U;                       //!< gauge configuration.

    std::vector<double> m_boundary_each_node; //!< b.c. on each node.

    //! arrays for communication buffer.
    double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
    double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
    double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
    double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

    Field m_w1, m_w2;   //!< working fields

   public:
    //! standard constructor.
    Fopr_WilsonGeneral(const Parameters& params) { init(params); }

    DEPRECATED
    Fopr_WilsonGeneral() { init("Dirac"); }

    DEPRECATED
    Fopr_WilsonGeneral(const std::string repr) { init(repr); }

    ~Fopr_WilsonGeneral() { tidyup(); }

    void set_parameters(const Parameters& params);

    void set_parameters(const double kappa_s, const double kappa_t,
                        const double nu_s, const double r_s,
                        const std::vector<int> bc);

    void get_parameters(Parameters& params) const;

    void set_config(Field *U);

    void set_mode(const std::string mode);

    std::string get_mode() const { return m_mode; }

    void mult(Field& v, const Field& w);

    void mult_dag(Field& v, const Field& w);

    void mult(Field& v, const Field& w, const std::string mod);

    void mult_dag(Field& v, const Field& w, const std::string mode);

    void mult_gm5(Field& v, const Field& w);

    void D(Field& v, const Field& w);

    void Ddag(Field& v, const Field& w);

    void DdagD(Field& v, const Field& w);

    void DDdag(Field& v, const Field& w);

    void H(Field& v, const Field& w);

    void mult_up(const int mu, Field&, const Field&);

    void mult_dn(const int mu, Field&, const Field&);

    int field_nvol() { return m_Nvol; }
    int field_nin()  { return m_Nvc * m_Nd; }
    int field_nex()  { return 1; }

    double flop_count();

   private:
    // prohibit copy
    Fopr_WilsonGeneral(const Fopr_WilsonGeneral&) {}
    Fopr_WilsonGeneral& operator=(const Fopr_WilsonGeneral&);

    //! obsolete initial setup.
    void init(const std::string repr);

    //! standard initial setup.
    void init(const Parameters& params);

    //! initial setup main.
    void setup();

    //! final clean-up.
    void tidyup();

    void mult_gm5_chiral(Field&, const Field&);
    void mult_gm5_dirac(Field&, const Field&);

    void D_chiral(Field&, const int ex1, const Field&, const int ex2);
    void D_dirac(Field&, const int ex1, const Field&, const int ex2);

    void mult_xp(Field&, const Field&);
    void mult_xm(Field&, const Field&);
    void mult_yp(Field&, const Field&);
    void mult_ym(Field&, const Field&);
    void mult_zp(Field&, const Field&);
    void mult_zm(Field&, const Field&);

    void mult_tp_dirac(Field&, const Field&);
    void mult_tm_dirac(Field&, const Field&);
    void mult_tp_chiral(Field&, const Field&);
    void mult_tm_chiral(Field&, const Field&);

    void daxpy(Field&, const double, const Field&);
    void clear(Field&);


#ifdef USE_FACTORY
   private:
    static Fopr *create_object() { return new Fopr_WilsonGeneral(); }

    static Fopr *create_object_with_repr(const std::string& repr)
    { return new Fopr_WilsonGeneral(repr); }

    static Fopr *create_object_with_params(const Parameters& params)
    { return new Fopr_WilsonGeneral(params); }

   public:
    static bool register_factory()
    {
      bool init = true;
      init &= Fopr::Factory_noarg::Register("WilsonGeneral/Imp", create_object);
      init &= Fopr::Factory_string::Register("WilsonFeneral/Imp",
                                             create_object_with_repr);
      init &= Fopr::Factory_params::Register("WilsonGeneral/Imp",
                                             create_object_with_params);
      return init;
    }
#endif
  };
}
#endif /* FOPR_WILSON_GENERAL_IMPL_IMP_INCLUDED */
