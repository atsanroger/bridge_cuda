/*!
        @file    fopr_Clover_eo.h

        @brief

        @author  Satoru Ueda (maintailed by H. Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-11-13 13:02:59 #$

        @version $LastChangedRevision: 2553 $
*/

#ifndef FOPR_CLOVER_EO_INCLUDED
#define FOPR_CLOVER_EO_INCLUDED

#include "Fopr/fopr_Wilson_eo.h"
#include "Fopr/fopr_CloverTerm_eo.h"
#include "Field/field_F.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Even-odd Clover fermion operator.

/*!
    This class is an even-odd version of Clover fermion operator.
    At present this is rough implementation, while correctly
    works, and to be updated by supplying complete functionality.
    Only the functions needed for even-odd preconditioned solver
    are ready.
                                        [20 June 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    Modify this code to work.           [03 Mar 2013 Y.Namekawa]
    Multi-threaded.                     [12 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
    Modifed toward ver.2.0.
                                        [04 Dec 2021 H.Matsufuru]
 */

class Fopr_Clover_eo : public Fopr_eo
{
 public:
  static const std::string class_name;

 private:
  // input parameters
  double m_kappa;               //!< hopping parameter.
  double m_cSW;                 //!< clover coefficient.
  std::vector<int> m_boundary;  //!< boundary condition.
  std::string m_repr;           //!< gamma-matrix type
  Bridge::VerboseLevel m_vl;    //!< verbose level

  std::string m_mode;

  // internal data members
  int m_Nc, m_Nd, m_NinF, m_Ndim;
  int m_Nvol, m_Nvol2;

  Fopr_Wilson_eo *m_fopr_w;
  Fopr_CloverTerm_eo *m_fopr_csw;

  Index_eo m_idx;
  Field_G *m_U;

  Field m_w1, m_w2;       //!< working field
  Field m_v1;             //!< working field

  bool m_is_initial_step; //!< to avoid redundant setup

 public:
  DEPRECATED
  Fopr_Clover_eo(const std::string repr) { init(repr); }

  Fopr_Clover_eo(const Parameters& params) { init(params); }

  ~Fopr_Clover_eo() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f);

  //- method for even odd fermion operator
  void preProp(Field& Be, Field& bo, const Field& b);

  void postProp(Field& x, const Field& xe, const Field& bo);

  //  const Field_F mult_csw_inv(const Field_F&, const int ieo);

  void mult_gm5(Field&, const Field&);

  void D(Field& v, const Field& f);
  void Ddag(Field& v, const Field& f);
  void DdagD(Field& v, const Field& f);
  void DDdag(Field& v, const Field& f);
  void H(Field& v, const Field& f);

  // Meo: ieo = 0: even <- odd, ieo=1: odd <- even
  void Meo(Field&, const Field&, const int ieo);
  void Mdageo(Field&, const Field&, const int ieo);

  void mult_isigma(Field_F& w, const Field_F& f,
                   const int mu, const int nu);

  int field_nin()  { return 2 * m_Nc * m_Nd; }
  int field_nvol() { return m_Nvol2; }
  int field_nex()  { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init(const std::string repr);

  void init(const Parameters& params);

  void setup();

  void tidyup();

#ifdef USE_FACTORY
 private:
  static Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Clover_eo(repr);
  }

  static Fopr *create_object_with_params(const Parameters& params)
  {
    return new Fopr_Clover_eo(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_string::Register("Clover_eo",
                                           create_object_with_repr);
    init &= Fopr::Factory_params::Register("Clover_eo",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
