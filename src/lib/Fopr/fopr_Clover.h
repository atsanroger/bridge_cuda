/*!
        @file    fopr_Clover.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_CLOVER_INCLUDED
#define FOPR_CLOVER_INCLUDED

#include "fopr_Wilson.h"
#include "fopr_CloverTerm.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover fermion operator.

/*!
    This class implements the clover (improved Wilson) fermion
    operator.
    Wilson kernel and clover term are implemented in other classes,
    and this class holds them as objects.
    (The implementation was modified after revision 645: before
     that, clover term was being implemented inside this class.)
    The `mode' controls which of D, Ddag, H, DdagD are multiplied
    when mult or mult_dag is called.
       [first ver. 24 Dec 2011/ modified 28 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                 [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                         [21 Mar 2015 Y.Namekawa]
    Modified toward ver.2.0.
                                        [04 Dec 2021 H.Matsufuru]
 */

class Fopr_Clover : public Fopr
{
 public:
  static const std::string class_name;

 private:
  double m_kappa;               //!< hopping parameter
  double m_cSW;                 //!< clover coefficient
  std::vector<int> m_boundary;  //!< boundary conditions
  std::string m_repr;           //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;    //!< verbose level

  std::string m_mode;           //!<  mode of multiplication

  // internal data members
  int m_Nc, m_Nd, m_NinF;       //!< internal parameters
  int m_Nvol, m_Ndim;           //!< internal parameters

  Fopr_Wilson *m_fopr_w;        //!< Wilson fermion kernel
  Fopr_CloverTerm *m_fopr_csw;  //!< Clover term operator
  const Field_G *m_U;           //!< gauge configuration (pointer)

  Field m_v1, m_v2;             //!< working field.

  bool m_is_initial_step;       //!< to avoid redundant setup

 public:
  DEPRECATED
  Fopr_Clover() { init("Dirac"); }

  DEPRECATED
  Fopr_Clover(const std::string repr) { init(repr); }

  Fopr_Clover(const Parameters& params) { init(params); }

  ~Fopr_Clover() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const double cSW,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f);

  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void DdagD(Field&, const Field&);
  void DDdag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w);

  void mult_up(const int mu, Field& v, const Field& w);

  void mult_dn(const int mu, Field& v, const Field& w);

  void mult_isigma(Field_F&, const Field_F&,
                   const int mu, const int nu);

  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:

  void init(const std::string repr);

  void init(const Parameters& params);

  void setup();

  void tidyup();


#ifdef USE_FACTORY
 private:
  static Fopr *create_object() { return new Fopr_Clover(); }

  static Fopr *create_object_with_arg(const std::string& repr)
  { return new Fopr_Clover(repr); }

  static Fopr *create_object_with_params(const Parameters& params)
  { return new Fopr_Clover(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Clover", create_object);
    init &= Fopr::Factory_string::Register("Clover",
                                           create_object_with_arg);
    init &= Fopr::Factory_params::Register("Clover",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
