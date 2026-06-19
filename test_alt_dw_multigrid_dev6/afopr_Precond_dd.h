/*!
      @file    afopr_Precond_dd.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2604 $
*/

#ifndef AFOPR_PRECOND_DD_INCLUDED
#define AFOPR_PRECOND_DD_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Parameters/commonParameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt/Solver/asolver.h"

class Field;

template<typename AFIELD>
class AFopr_Precond_dd : public AFopr_dd<AFIELD>
{
 public:
  typedef typename AFIELD::real_t      real_t;
  typedef typename AFIELD::complex_t   complex_t;
  static const std::string class_name;

 protected:

  AFopr_dd<AFIELD> *m_afopr;
  ASolver<AFIELD> *m_asolver;

  AFIELD m_v1, m_w1, m_w2;

  std::string m_mode;   //!< mult mode
  Bridge::VerboseLevel m_vl;  //!< verbose level


 public:
  //! constructor.
  AFopr_Precond_dd(AFopr_dd<AFIELD>* fopr,
		   ASolver<AFIELD>* solver,
		   const Parameters& params)
  { init(fopr, solver, params); }

  //! destructor.
  ~AFopr_Precond_dd() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);


  //! setting gauge configuration.
  void set_config(Field *u);

  //! QXS version requires convert of spinor field.
  bool needs_convert() { return m_afopr->needs_convert(); }

  //! convert of spinor field.
  void convert(AFIELD& v, const Field& w){
    m_afopr->convert(v,w);
  }

  //! reverse of spinor field.
  void reverse(Field& v, const AFIELD& w){
    m_afopr->reverse(v,w);
  }

  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult(AFIELD&, const AFIELD&, const string mode);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);


  void mult_dup(AFIELD&, const AFIELD&, int mu);
  void mult_ddn(AFIELD&, const AFIELD&, int mu);

  void mult_sap(AFIELD&, const AFIELD&, const int eo);

  void mult_dd(AFIELD&, const AFIELD&);


  //! returns inner size parameter.
  int field_nin() { return m_afopr->field_nin(); }

  //! returns local volume size parameter.
  int field_nvol() { return m_afopr->field_nvol(); }

  //! returns external size parameter.
  int field_nex() { return m_afopr->field_nex(); }

  //! returns floating operation counts.
  double flop_count() { return flop_count(m_mode); }

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

  double flop_count_sap();

  // for debug
  void dump();

 private:
  //! initial setup.
  void init(AFopr_dd<AFIELD>* fopr, ASolver<AFIELD>* solver,
            const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  //  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  void mult_D_sap(AFIELD&, const AFIELD&, const int eo);
  void mult_Ddag_sap(AFIELD&, const AFIELD&, const int eo);

  void mult_gm4(AFIELD&, const AFIELD&);

  //! standard D mult.
  void mult_D(AFIELD&, const AFIELD&);

  void clear(real_t *);

  //  void aypx(complex_t, real_t *, real_t *);
  //  void gm5_aypx(complex_t, real_t *, real_t *);


#ifdef USE_FACTORY
  // For this class, temporarily factory is not defined.
  /*
 private:
  static AFopr<AFIELD> *create_object_with_fopr_solver_params(
                                        AFopr_dd<AFIELD>* fopr,
                                        ASolver<AFIELD>* solver,
                                        const Parameters& params)
  { return new AFopr_Precond_dd(fopr, solve, params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
            "Precond_dd", create_object_with_fopr_solver_params);
    return init1;
  }
  */
#endif
};

#endif  // AFOPR_PRECOND_DD_INCLUDED
