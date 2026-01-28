/*!
        @file    fopr_Wilson_eo_impl.h

        @brief

        @author  Satoru UEDA (maintailed by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_EO_IMPL_IMP_INCLUDED
#define FOPR_WILSON_EO_IMPL_IMP_INCLUDED

#include "Fopr/fopr_eo.h"
#include "Field/index_eo.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Implementation of even-odd Wilson fermion operator.

/*!
    This class is a subclass of Fopr_Wilson_eo and implements
    an even-odd version of the Wilson fermion operator.
    This is rather straightforward and readable while slower
    version that was coded by S.Ueda [20 Jun 2012 S.UEDA].
    The implementation class was separated.
                                     [07 Jul 2014 H.Matsufuru]
    Multi-threading method is changed and applied to methods
    other than mult.
                                     [03 Jul 2021 H.Matsufuru]
 */

namespace Imp {
  class Fopr_Wilson_eo : public Fopr_eo
  {
   public:
    static const std::string class_name;

   private:
    // input parameters
    double m_kappa;                //!< hopping parameter
    std::vector<int> m_boundary;   //!< boundary condition
    std::string m_repr;            //!< Dirac matrix representation
    Bridge::VerboseLevel m_vl;     //!< verbose level

    std::string m_mode;            //!< mult mode.

    // internal data members
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nvol, m_Nvol2, m_Ndim;
    int m_Nx, m_Ny, m_Nz, m_Nt, m_Nx2;

    std::vector<double> m_boundary_each_node; //!< b.c. for each node.

    Index_eo m_index;
    Field_G m_Ueo;

    std::vector<int> m_yzt_eo;   //!< yzt parity
    Field m_v1, m_v2;            //!< working field
    Field m_w1, m_w2;            //!< working field

    //! communication buffers
    double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
    double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
    double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
    double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

   public:
    DEPRECATED
    Fopr_Wilson_eo() { init("Dirac"); }

    DEPRECATED
    Fopr_Wilson_eo(const std::string repr) { init(repr); }

    Fopr_Wilson_eo(const Parameters& params) { init(params); }

    ~Fopr_Wilson_eo() { tidyup(); }

    void set_parameters(const Parameters& params);

    void set_parameters(const double kappa, const std::vector<int> bc);

    void get_parameters(Parameters& params) const;

    void set_config(Field *U);
    void set_config_omp(Field *U);
    void set_config_impl(Field *U);

    void set_mode(const std::string mode);

    std::string get_mode() const { return m_mode; }

    // method for even odd fermion operator
    void preProp(Field& Be, Field& bo, const Field& b);

    void postProp(Field& x, const Field& xe, const Field& bo);

    void mult(Field& v, const Field& w);

    void mult_dag(Field& v, const Field& w);

    void mult(Field& v, const Field& w, const std::string mode);

    void mult_dag(Field& v, const Field& w, const std::string mode);

    void D(Field& v, const Field& w);
    void Ddag(Field& v, const Field& w);
    void DdagD(Field& v, const Field& w);
    void DDdag(Field& v, const Field& w);
    void H(Field& v, const Field& w);

    // Meo: ieo = 0: even <- odd, 1: odd <- even
    void Meo(Field&, const Field&, const int ieo);

    void Mdageo(Field&, const Field&, const int ieo);

    //  void MeoMoe(Field&, const Field&);
    void Meo_gm5(Field&, const Field&, const int ieo);

    void mult_gm5(Field&, const Field&);
    void mult_gm5(Field&);
    void gm5_dirac(Field&, const Field&);
    void gm5_chiral(Field&, const Field&);

    //! gamma_5 (1 - gamma_mu) v(x + mu) used in force calculation.
    void gm5p(const int mu, Field&, const Field& v);

    int field_nvol() { return m_Nvol2; }
    int field_nin()  { return m_Nvc * m_Nd; }
    int field_nex()  { return 1; }

    //! this returns the number of floating point operations of Meo.
    double flop_count();

   private:
    void init(const std::string);

    void init(const Parameters&);

    void setup();

    void tidyup();

    void Meo_dirac(Field&, const Field&, const int ieo);
    void Meo_chiral(Field&, const Field&, const int ieo);

    void mult_xp(Field&, const Field&, const int ieo);
    void mult_xm(Field&, const Field&, const int ieo);
    void mult_yp(Field&, const Field&, const int ieo);
    void mult_ym(Field&, const Field&, const int ieo);
    void mult_zp(Field&, const Field&, const int ieo);
    void mult_zm(Field&, const Field&, const int ieo);
    void mult_tp_dirac(Field&, const Field&, const int ieo);
    void mult_tm_dirac(Field&, const Field&, const int ieo);
    void mult_tp_chiral(Field&, const Field&, const int ieo);
    void mult_tm_chiral(Field&, const Field&, const int ieo);


#ifdef USE_FACTORY
   private:
    static Fopr *create_object()
    {
      return new Fopr_Wilson_eo();
    }

    static Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson_eo(repr);
    }

    static Fopr *create_object_with_params(const Parameters& params)
    {
      return new Fopr_Wilson_eo(params);
    }

   public:
    static bool register_factory()
    {
      bool init = true;
      init &= Fopr::Factory_noarg::Register("Wilson_eo/Imp",
                                            create_object);
      init &= Fopr::Factory_string::Register("Wilson_eo/Imp",
                                             create_object_with_repr);
      init &= Fopr::Factory_params::Register("Wilson_eo/Imp",
                                             create_object_with_params);
      return init;
    }
#endif
  };
}
#endif /* FOPR_WILSON_EO_IMPL_IMP_INCLUDED */
