/*!
        @file    afopr_Overlap.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$

        @version $LastChangedRevision: 2530 $
*/


#ifndef AFOPR_OVERLAP_INCLUDED
#define AFOPR_OVERLAP_INCLUDED

#include "Fopr/afopr_Sign.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Overlap fermion operator.

/*!
    This class defines the overlap fermion operator with a given
    kernel fermion operator that currently is assumed as the Wilson
    fermion operator.
    All the necessary parameters for the kernel fermion operator must
    be set outside of the class, while the parameters for the sign
    function operator are set in this class.
    When eigenvalues and eigenvectors of the base fermion operator
    are given, these modes are used in the sign function to improve
    the approximation formula.
                                     [20 Dec 2011 H.Matsufuru]
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
    Changed to a template class in ver.2.0.
                                     [12 Feb 2022 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Overlap : public AFopr<AFIELD>
{
 public:
  typedef AFopr<AFIELD>             AFOPR;
  typedef typename AFIELD::real_t   real_t;
  static const std::string class_name;

 private:

  // input parameters
  real_t m_mq;                 //!< quark mass
  real_t m_M0;                 //!< domain-wall height
  int m_Np;                    //!< number of poles in rational approx.
  real_t m_x_min;              //!< valid range of approximate sign function
  real_t m_x_max;              //!< valid range of approximate sign function
  int m_Niter;                 //!< max iteration of shiftsolver
  double m_Stop_cond;          //!< stopping condition of shift solver
  std::vector<int> m_boundary; //!< boundary condition

  std::string m_kernel_type;   //!< kernel type (if given)
  std::string m_repr;          //!< gamma-matrix repr. (if given)

  Bridge::VerboseLevel m_vl;

  std::string m_mode;  //!< mult mode

  // local variables
  int m_Nin, m_Nvol, m_Nex;

  AFOPR *m_fopr_w;            //!< kernel fermion operator.
  bool m_kernel_created;      //!< whether kernel is created in this object

  AFopr_Sign<AFIELD> *m_sign; //!< sign function approximation.

  //- for low-eigenmode subtraction.
  int m_Nsbt;
  std::vector<real_t> *m_ev;
  std::vector<AFIELD> *m_vk;

  //- workarea
  AFIELD m_w0, m_w1, m_w2;

  bool m_is_initial_step; //!< to avoid redundant setup

 public:

  AFopr_Overlap(const Parameters& params) { init(params); }

  AFopr_Overlap(AFOPR *fopr, const Parameters& params)
  { init(fopr, params); }

  DEPRECATED
  AFopr_Overlap(AFOPR *fopr)
    : m_fopr_w(fopr) { init(); }

  ~AFopr_Overlap() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const real_t mq, const real_t M0, const int Np,
                      const real_t x_min, const real_t x_max,
                      const int Niter, const real_t Stop_cond,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_lowmodes(const int Nsbt, std::vector<real_t> *,
                    std::vector<AFIELD> *);

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  void H(AFIELD& v, const AFIELD& w);
  void D(AFIELD& v, const AFIELD& w);
  void Ddag(AFIELD& v, const AFIELD& w);
  void DdagD(AFIELD& v, const AFIELD& w);

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr_w->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr_w->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr_w->reverse(v, w); }

  int field_nin() { return m_Nin; }
  int field_nvol() { return m_Nvol; }
  int field_nex() { return m_Nex; }

  //! returns the number of floating point operations.
  double flop_count();

  //! returns the flops per site for specified mode.
  double flop_count(const std::string mode);

 private:

  void init(const Parameters& params);

  void init(AFOPR *fopr, const Parameters& params);

  void init();

  void tidyup();


#ifdef USE_FACTORY
 private:
  static Fopr *create_object(AFOPR *fopr)
  { return new AFopr_Overlap<AFIELD>(fopr); }

  static Fopr *create_object_params(const Parameters& params)
  { return new AFopr_Overlap<AFIELD>(params); }

  static Fopr *create_object_fopr_params(AFOPR *fopr,
                                         const Parameters& params)
  { return new AFopr_Overlap<AFIELD>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AFOPR::Factory_fopr::Register("Overlap", create_object);
    init &= AFOPR::Factory_params::Register("Overlap",
                                            create_object_params);
    init &= AFOPR::Factory_fopr_params::Register("Overlap",
                                                 create_object_fopr_params);
    return init;
  }
#endif
};
#endif
