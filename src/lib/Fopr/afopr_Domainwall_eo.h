/*!
        @file    afopr_Domainwall_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-06-08 10:27:09 #$

        @version $LastChangedRevision: 2522 $
*/

#ifndef AFOPR_DOMAINWALL_EO_INCLUDED
#define AFOPR_DOMAINWALL_EO_INCLUDED

#include <vector>
#include <string>

#include "lib/Fopr/afopr_eo.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

class Field;
class Field_G;
template<typename AFIELD>
class Index_eo_Domainwall;

//! Domain-wall fermion operator with even-odd site index.

/*!
   This class implements the even-odd version of Domain-wall
   fermion including the Mobius form.
                                        [18 Apr 2017 H.Matsufuru]
   This class was developed as a template class in the alternative
   branch and then incorporated into ver.2.0.
                                        [07 Mar 2022 H.Matsufuru]
   H.Neff's improving parameter alpha is incorporated.
                                        [13 Jan 2025 H.Matsufuru]
 */
template<typename AFIELD>
class AFopr_Domainwall_eo : public AFopr_eo<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                 //!< quark mass
  real_t m_M0;                 //!< domain-wall height
  int m_Ns;                    //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<real_t> m_b;
  std::vector<real_t> m_c;
  real_t m_alpha;
  Bridge::VerboseLevel m_vl; //!< verbose level
  std::string m_kernel_type;
  std::string m_repr;

  std::string m_mode;

  int m_NinF;                //!< on-site d.o.f.
  int m_Nvol;                //!< volume size of even or odd vector.
  int m_Nvol2;               //!< volume size of even or odd vector.
  int m_Ndim;                //!< spacetime dimensions

  Index_eo_Domainwall<AFIELD> *m_index_eo;

  AFopr<AFIELD> *m_foprw;

  AFIELD m_w1, m_v1, m_v2;        //!< woking 5d vectors
  AFIELD m_w4, m_v4, m_y4, m_t4;  //!< woking 4d vectors

  // for convert and reverse
  Field m_w4lex;
  AFIELD m_v4lex;


  // for preconditioning
  std::vector<real_t> m_dp;
  std::vector<real_t> m_dm;
  std::vector<real_t> m_e;
  std::vector<real_t> m_f;
  real_t m_g;

  //! initial setup.
  void init(const Parameters& params);

  //! final tidyup.
  void tidyup();

 public:
  //! constructor.
  AFopr_Domainwall_eo(const Parameters& params)
    : AFopr_eo<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Domainwall_eo() { tidyup(); }

  void set_parameters(const Parameters& params);

  //! set parameters in the case of Moebius domain-wall.
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns, const std::vector<int> bc,
                      const real_t b, const real_t c,
                      const real_t alpha);

  void get_parameters(Parameters& params) const;

  //! set parameters of kernel operaotr.
  void set_kernel_parameters(const Parameters& params);

  //! set parameters for preconditioning.
  void set_precond_parameters();

  //! set coefficients if they depend in s.
  void set_coefficients(const std::vector<real_t> b,
                        const std::vector<real_t> c);

  //! this class needs convert of fermion field.
  bool needs_convert() { return m_foprw->needs_convert(); }

  //! convert Field to AField for this class.
  void convert(AFIELD&, const Field&);

  //! reverse AField to Field.
  void reverse(Field&, const AFIELD&);

  void set_config(Field *U) { m_foprw->set_config(U); }

  void set_config(unique_ptr<Field_G>& U)
  { m_foprw->set_config(U.get()); }

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  // methods for even odd fermion operator
  void preProp(AFIELD& Be, AFIELD& bo, const AFIELD& b);

  void postProp(AFIELD& x, const AFIELD& xe, const AFIELD& bo);

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w);

  void mult(AFIELD& v, const AFIELD& w, const std::string mode);

  void mult_gm5(AFIELD& v, const AFIELD& w);

  void mult_gm5_4d(AFIELD& v, const AFIELD& w);


  void DdagD(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void Hdag(AFIELD&, const AFIELD&);

  void D_ee(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_ee(AFIELD&, const AFIELD&, const int ieo);
  void D_eo(AFIELD&, const AFIELD&, const int ieo);
  void Ddag_eo(AFIELD&, const AFIELD&, const int ieo);

  void mult_gm5R(AFIELD&, const AFIELD&);
  void mult_R(AFIELD&, const AFIELD&);

  //  preconditioner
  void L_inv(AFIELD&, const AFIELD&);
  void U_inv(AFIELD&, const AFIELD&);
  void Ldag_inv(AFIELD&, const AFIELD&);
  void Udag_inv(AFIELD&, const AFIELD&);

  int field_nin() { return m_NinF; }
  int field_nvol() { return m_Nvol2; }
  int field_nex() { return m_Ns; }

  //! this returns the number of floating point number operations.
  double flop_count() { return flop_count(m_mode); }

  //! flop-count for specified mode.
  double flop_count(std::string mode);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Domainwall_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
      "Domainwall_eo", create_object_with_params);
    return init1;
  }
#endif
};

#endif
