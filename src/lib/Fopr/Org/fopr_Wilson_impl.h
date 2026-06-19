/*!
        @file    fopr_Wilson_impl.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_WILSON_IMPL_ORG_INCLUDED
#define FOPR_WILSON_IMPL_ORG_INCLUDED

#include "Fopr/fopr.h"

#include "Field/field_F.h"
#include "Field/shiftField_lex.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator.

/*!
    This class implements the standard Wilson fermion operator
    by making use of general tools available in this code set.
    Original version was written in 24 Dec 2011.
    The implemetation was modified to fit the up-to date style
    as well as to enable multi-threading.
                                    [17 Feb 2022 H.Matsufuru]
 */

namespace Org {
  class Fopr_Wilson : public Fopr
  {
   public:
    static const std::string class_name;

   private:

    // input parameters
    double m_kappa;              //!< hopping parameter
    std::vector<int> m_boundary; //!< boundary condition
    std::string m_repr;          //!< Dirac matrix representation
    Bridge::VerboseLevel m_vl;   //!< verbose level

    std::string m_mode;          //!< matrix multiplcation mode

    // local variables
    int m_Nc, m_Nd, m_Nvol, m_Ndim;

    const Field_G *m_U;     //!< pointer ot gauge configuration

    std::vector<GammaMatrix> m_GM;

    ShiftField_lex *m_shift;

    Field m_w1, m_w2, m_w3;   //!< working fields
    Field_F m_v1, m_v2;       //!< working fields with spinor structure

   public:

    //! standard constructor with input parameters.
    Fopr_Wilson(const Parameters& params) { init(params); }

    DEPRECATED
    Fopr_Wilson() { init("Dirac"); }

    DEPRECATED
    Fopr_Wilson(const std::string repr) { init(repr); }

    ~Fopr_Wilson() { tidyup(); }

    void set_parameters(const Parameters& params);

    void set_parameters(const double kappa, const std::vector<int> bc);

    void get_parameters(Parameters& params) const;

    void set_config(Field *U);

    void set_mode(const std::string mode);

    std::string get_mode() const { return m_mode; }

    void mult(Field& v, const Field& w);

    void mult_dag(Field& v, const Field& w);

    void mult(Field& v, const Field& w, const std::string mod);

    void mult_dag(Field& v, const Field& w, const std::string mode);

    void mult_gm5(Field& v, const Field& w);

    void mult_up(const int mu, Field& v, const Field& w);

    void mult_dn(const int mu, Field& v, const Field& w);

    void D(Field&, const Field&);

    void D_ex(Field& v, const int ex1, const Field& w, const int ex2);

    void Ddag(Field&, const Field&);

    void DdagD(Field&, const Field&);

    void DDdag(Field&, const Field&);

    void H(Field&, const Field&);

    void proj_chiral(Field& w, const int ex1, const Field& v,
                     const int ex2, const int ipm);

    void mult_gm5p(const int mu, Field& v, const Field& w);

    int field_nin()  { return 2 * m_Nc * m_Nd; }
    int field_nvol() { return m_Nvol; }
    int field_nex()  { return 1; }

    //! this returns the number of floating point operations.
    double flop_count();

   private:
    //- prohibit copy
    Fopr_Wilson(const Fopr_Wilson&) {}
    Fopr_Wilson& operator=(const Fopr_Wilson&);

    //! standard initial setup.
    void init(const Parameters& params);

    //! obsolete initial setup.
    void init(std::string repr);

    //! common parts of setup.
    void setup();

    //! final clean-up.
    void tidyup();


#ifdef USE_FACTORY
   private:
    static Fopr *create_object()
    {
      return new Fopr_Wilson();
    }

    static Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson(repr);
    }

    static Fopr *create_object_with_params(const Parameters& params)
    {
      return new Fopr_Wilson(params);
    }

   public:
    static bool register_factory()
    {
      bool init = true;
      init &= Fopr::Factory_noarg::Register("Wilson/Org", create_object);
      init &= Fopr::Factory_string::Register("Wilson/Org", create_object_with_repr);
      init &= Fopr::Factory_params::Register("Wilson/Org", create_object_with_params);
      return init;
    }
#endif
  };
}
#endif /* FOPR_WILSON_IMPL_ORG_INCLUDED */
