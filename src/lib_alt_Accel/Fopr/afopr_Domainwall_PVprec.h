/*!
        @file    afopr_Domainwall_PVprec.h

        @brief   Pauli-Villars + site-diagonal preconditioned Moebius
                 domain-wall operator.

        @author  Wei-Lun Chen

        @details
          Implements the preconditioned operator

              A = (D_PV C_PV^{-1})^dagger  C^{-1} D

          which was identified as the best-conditioned target operator
          for multigrid in
            I. Kanamori, W.-L. Chen, H. Matsufuru,
            "Spectrum of Preconditioned Moebius Domain-wall Operators",
            PoS(LATTICE2024)414.

          Here
            D       : Moebius domain-wall operator at the physical mass,
            C       : its site-diagonal (5th-dimension) part, so that
                      C^{-1} is the exact LU inverse exposed as mode "Prec",
            D_PV    : the Pauli-Villars operator, i.e. the same operator at
                      quark_mass = quark_mass_PauliVillars (default 1.0),
            C_PV    : the site-diagonal part of D_PV (mode "Prec" of the PV
                      instance).

          All four factors are applied exactly; no iterative solver is used.
          This is the FP (single/double precision) reference implementation.
*/

#ifndef AFOPR_DOMAINWALL_PVPREC_INCLUDED
#define AFOPR_DOMAINWALL_PVPREC_INCLUDED

#include "lib/Fopr/afopr.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

template<typename AFIELD>
class AFopr_Domainwall_PVprec : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  AFopr_Domainwall_5din<AFIELD> *m_fopr;     //!< physical-mass operator D
  AFopr_Domainwall_5din<AFIELD> *m_fopr_PV;  //!< Pauli-Villars operator D_PV

  std::string m_mode;                        //!< mult mode ("D","Ddag","DdagD")

  AFIELD m_t1, m_t2, m_t3;                    //!< working 5d vectors

 public:
  //! constructor.
  AFopr_Domainwall_PVprec(const Parameters& params) : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Domainwall_PVprec() { tidyup(); }

  void set_parameters(const Parameters&) {}

  //! set gauge configuration on both internal operators.
  void set_config(Field *U)
  {
    m_fopr->set_config(U);
    m_fopr_PV->set_config(U);
  }

  void set_mode(std::string mode) { m_mode = mode; }
  std::string get_mode() const { return m_mode; }

  void mult(AFIELD& v, const AFIELD& w)
  {
    if (m_mode == "DdagD") {
      D(m_t3, w);
      Ddag(v, m_t3);
    } else if (m_mode == "Ddag") {
      Ddag(v, w);
    } else {
      D(v, w);
    }
  }

  void mult_dag(AFIELD& v, const AFIELD& w)
  {
    if (m_mode == "DdagD") {
      D(m_t3, w);
      Ddag(v, m_t3);
    } else if (m_mode == "Ddag") {
      D(v, w);
    } else {
      Ddag(v, w);
    }
  }

  void mult(AFIELD& v, const AFIELD& w, std::string mode)
  {
    if (mode == "DdagD") {
      D(m_t3, w);
      Ddag(v, m_t3);
    } else if (mode == "Ddag") {
      Ddag(v, w);
    } else if (mode == "leftprec") {
      mult_leftprec(v, w);     // v = P_L w = (D_PV C_PV^{-1})^dag C^{-1} w
    } else if (mode == "physicalD") {
      mult_physical_D(v, w);   // v = D w (physical Moebius, for residual check)
    } else {
      D(v, w);
    }
  }

  //! v = P_L w = (D_PV C_PV^{-1})^dagger C^{-1} w  (the left preconditioner).
  //! Use to build the RHS  b' = P_L eta  so that solving  A psi = b'  yields
  //! psi = D^{-1} eta exactly (A = P_L D), i.e. the SAME propagator as a direct
  //! D psi = eta solve (CGNE). Forms b' with one application; uses m_t1/m_t2,
  //! which are free outside an A-mult, so this is safe to call before solving.
  void mult_leftprec(AFIELD& v, const AFIELD& w)
  {
    m_fopr->mult(m_t2, w, "Prec");        // t2 = C^{-1} w
    m_fopr_PV->mult_dag(m_t1, m_t2, "D"); // t1 = D_PV^dagger t2
    m_fopr_PV->mult_dag(v, m_t1, "Prec"); // v  = C_PV^{-dagger} t1
  }

  //! v = D w : the physical-mass Moebius operator only (no preconditioning).
  //! For the true-residual check ||D psi - eta|| that certifies the propagator.
  void mult_physical_D(AFIELD& v, const AFIELD& w)
  {
    m_fopr->mult(v, w, "D");
  }

  bool needs_convert() { return m_fopr->needs_convert(); }
  void convert(AFIELD& v, const Field& w) { m_fopr->convert(v, w); }
  void reverse(Field& v, const AFIELD& w) { m_fopr->reverse(v, w); }

  int field_nin()  { return m_fopr->field_nin(); }
  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nex()  { return m_fopr->field_nex(); }

  double flop_count() { return 0.0; }

 private:
  void init(const Parameters& params)
  {
    m_vl   = CommonParameters::Vlevel();
    m_mode = "D";

    // physical-mass operator D
    m_fopr = new AFopr_Domainwall_5din<AFIELD>(params);

    // Pauli-Villars operator D_PV: same parameters but quark_mass = M_PV.
    Parameters params_PV = params;
    double M_PV = 1.0;
    if (params.is_set("quark_mass_PauliVillars")) {
      M_PV = params.get_double("quark_mass_PauliVillars");
    }
    params_PV.set_double("quark_mass", M_PV);
    m_fopr_PV = new AFopr_Domainwall_5din<AFIELD>(params_PV);

    vout.general(m_vl, "%s: constructed,"
                 " Pauli-Villars mass M_PV = %f\n",
                 class_name.c_str(), M_PV);

    const int nin  = m_fopr->field_nin();
    const int nvol = m_fopr->field_nvol();
    const int nex  = m_fopr->field_nex();
    m_t1.reset(nin, nvol, nex);
    m_t2.reset(nin, nvol, nex);
    m_t3.reset(nin, nvol, nex);
  }

  void tidyup()
  {
    delete m_fopr_PV;
    delete m_fopr;
  }

  //! v = (D_PV C_PV^{-1})^dagger C^{-1} D w
  void D(AFIELD& v, const AFIELD& w)
  {
    m_fopr->mult(m_t1, w, "D");             // t1 = D w
    m_fopr->mult(m_t2, m_t1, "Prec");       // t2 = C^{-1} D w
    m_fopr_PV->mult_dag(m_t1, m_t2, "D");   // t1 = D_PV^dagger t2
    m_fopr_PV->mult_dag(v, m_t1, "Prec");   // v  = C_PV^{-dagger} t1
  }

  //! v = D^dagger C^{-dagger} (D_PV C_PV^{-1}) w   (adjoint of D)
  void Ddag(AFIELD& v, const AFIELD& w)
  {
    m_fopr_PV->mult(m_t1, w, "Prec");       // t1 = C_PV^{-1} w
    m_fopr_PV->mult(m_t2, m_t1, "D");       // t2 = D_PV t1
    m_fopr->mult_dag(m_t1, m_t2, "Prec");   // t1 = C^{-dagger} t2
    m_fopr->mult_dag(v, m_t1, "D");         // v  = D^dagger t1
  }

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_PVprec(params); }

 public:
  static bool register_factory()
  {
    return AFopr<AFIELD>::Factory_params::Register(
                     "Domainwall_PVprec", create_object_with_params);
  }
#endif
};

#endif // AFOPR_DOMAINWALL_PVPREC_INCLUDED
