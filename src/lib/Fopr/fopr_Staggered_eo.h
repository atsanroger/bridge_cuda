/*!
        @file    fopr_Staggered_eo.h
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_STAGGERED_EO_INCLUDED
#define FOPR_STAGGERED_EO_INCLUDED

#include "fopr_eo.h"

#include "bridge_defs.h"
#include "Field/index_eo.h"
#include "Field/shiftField_eo.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Even-odd staggered fermion operator.

/*!
    Standard staggered fermion operator with even-odd site index.
    For ver.2.0.                      [21 Nov 2021 H.Matsufuru]

    - Original version                [24 Dec 2011 H.Matsufuru]
    - Normalization of operator is changed from hopping to
      mass one.                                   [21 Nov 2021]
    - Largely modified toward ver.2.0             [21 Nov 2021]
 */

class Fopr_Staggered_eo : public Fopr_eo
{
 public:
  static const std::string class_name;

 private:
  // input parameters
  double m_mq;                 //!< fermion mass
  std::vector<int> m_boundary; //!< boundary conditions
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::string m_mode;          //!< mult mode

  // internal data members
  int m_Nc, m_Nvc;
  int m_Nvol, m_Nvol2, m_Ndim;

  Field_G m_Ueo;            //!< even-odd gauge configuration
  Field m_staggered_phase;  //!< staggered phase

  Index_eo m_index_eo;

  ShiftField_eo *m_shift_eo;

  Field m_t1, m_t2;    //!< working vectors
  Field m_v1, m_v2;    //!< working vectors

 public:

  DEPRECATED
  Fopr_Staggered_eo() { init(); }

  //! constructor with a Paramters object.
  Fopr_Staggered_eo(const Parameters& params) { init(params); }

  //! deconstructor.
  ~Fopr_Staggered_eo() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_config_impl(Field *U);

  void set_config_omp(Field *U);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(Field& w, const Field& f);

  void mult_dag(Field& w, const Field& f);

  void mult(Field& w, const Field& f, const std::string mode);

  void mult_dag(Field& w, const Field& f, const std::string mode);

  //! for even-odd solver: called before linear solver.
  void preProp(Field&, Field&, const Field&);

  //! for even-odd solver: called after linear solver.
  void postProp(Field&, const Field&, const Field&);

  //! D = m^2 - Deo * Doe.
  void MeoMoe(Field& v, const Field& f);

  //! hopping term: ieo = 0: even <-- odd, 1: odd  <-- even
  void Meo(Field&, const Field&, const int ieo, const int jd);

  //! multiplied staggered phase with a given vector.
  void mult_staggered_phase(Field&, const int mu, const int ieo);

  int field_nin() { return 2 * m_Nc; }

  int field_nvol() { return m_Nvol2; }

  int field_nex() { return 1; }

  //! returns the number of floating point operations.
  double flop_count() { return flop_count(m_mode); }

  //! returns the flop count for specified mode.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  void init();

  //! initial setup.
  void init(const Parameters& params);

  void setup();

  //! final clean-up.
  void tidyup();

  void set_staggered_phase();

  void mult_up(Field&, const Field&, const int mu, const int ieo);

  void mult_dn(Field&, const Field&, const int mu, const int ieo);


#ifdef USE_FACTORY
 private:
  static Fopr *create_object() { return new Fopr_Staggered_eo(); }

  static Fopr *create_object_with_params(const Parameters& params)
  { return new Fopr_Staggered_eo(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Fopr::Factory_noarg::Register("Staggered_eo", create_object);
    init &= Fopr::Factory_params::Register("Staggered_eo",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
