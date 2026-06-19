/*!
        @file    fopr_Wilson_eo_impl.h

        @brief

        @author  Satoru UEDA (maintailed by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_EO_IMPL_ORG_INCLUDED
#define FOPR_WILSON_EO_IMPL_ORG_INCLUDED

#include "Fopr/fopr_eo.h"
#include "Field/field_F.h"
#include "Field/shiftField_eo.h"
#include "Tools/gammaMatrixSet.h"
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
    Multi-threading is implemented.
    Implementation is updated suitably to ver.2.0.
                                     [20 Feb 2022 H.Matsufuru]
 */

namespace Org {
  class Fopr_Wilson_eo : public Fopr_eo
  {
   public:
    static const std::string class_name;

   private:
    //input parameters
    double m_kappa;               //!< hopping parameter
    std::vector<int> m_boundary;  //!< boundary condition
    std::string m_repr;           //!< gamma matrix representation
    Bridge::VerboseLevel m_vl;    //!< verbose level

    std::string m_mode;           //!< multiplication mode

    // local variables
    int m_Nvol, m_Nvol2, m_Ndim;
    int m_Nc, m_Nd;

    Field_G m_Ueo;                //!< even-odd configuration

    ShiftField_eo *m_shift;       //!< even-odd field shifter

    std::vector<GammaMatrix> m_GM;

    Index_eo m_index;

    Field_F m_v1, m_v2;  //!< working field (used in D, Ddag)
    Field_F m_v3;        //!< working field (used in DdagD, DDdag, H)

    Field_F m_t1, m_t2;  //!< working field (used in mult_up/dn)
    Field_F m_w1, m_w2;  //!< working field (used in Meo, Mdageo, mult_gm5)

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

    void set_mode(std::string mode);

    std::string get_mode() const { return m_mode; }

    void mult(Field&, const Field&);

    void mult_dag(Field&, const Field&);

    void mult(Field& v, const Field& w, const std::string mode);

    void mult_dag(Field& v, const Field& w, const std::string mode);

    void preProp(Field& Be, Field& bo, const Field& b);
    void postProp(Field& x, const Field& xe, const Field& bo);

    void D(Field&, const Field&);
    void Ddag(Field&, const Field&);
    void DdagD(Field&, const Field&);
    void DDdag(Field&, const Field&);
    void H(Field&, const Field&);

    // ieo=0: even <-- odd
    // ieo=1: odd  <-- even

    void Meo(Field&, const Field&, const int ieo);
    void Mdageo(Field&, const Field&, const int ieo);

    void mult_gm5(Field&, const Field&);

    //! gamma_5 (1 - gamma_mu) v(x + mu)
    void gm5p(const int mu, Field&, const Field& v);

    int field_nvol() { return m_Nvol2; }
    int field_nin() { return 2 * m_Nc * m_Nd; }
    int field_nex() { return 1; }

    //! this returns the number of floating point operations of Meo.
    double flop_count();

   private:
    //! initial setup (standard)
    void init(const Parameters&);

    //! initial setup (obsolete)
    void init(const std::string);

    //! common setup in initialization.
    void setup();

    //! final clean-up
    void tidyup();

    void mult_up(const int mu, Field_F&, const Field_F&, const int ieo);

    void mult_dn(const int mu, Field_F&, const Field_F&, const int ieo);

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
      init &= Fopr::Factory_noarg::Register("Wilson_eo/Org", create_object);
      init &= Fopr::Factory_string::Register("Wilson_eo/Org",
                                             create_object_with_repr);
      init &= Fopr::Factory_params::Register("Wilson_eo/Org",
                                             create_object_with_params);
      return init;
    }
#endif
  };
}
#endif /* FOPR_WILSON_EO_IMPL_ORG_INCLUDED */
