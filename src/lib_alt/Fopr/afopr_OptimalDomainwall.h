/*!                                                                             
        @file    afopr_OptimalDomainwall.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

#ifndef AFOPR_OPTIMALDOMAINWALL_INCLUDED
#define AFOPR_OPTIMALDOMAINWALL_INCLUDED

#include <valarray>
#include <string>

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Field/field.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Fopr/afopr_Domainwall_General.h"

class Field_G;

//! Alternative version of the Optimal Domain-wall fermion operator.
/*!
   This class implements optimal domain-wall fermion operator.
                                         [06 May 2017 H.Matsufuru]
 */
template<typename AFIELD>
class AFopr_OptimalDomainwall : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  double m_mq;            //!< quark mass
  double m_M0;            //!< domain-wall height
  int    m_Ns;            //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<double> m_b;
  std::vector<double> m_c;
  Bridge::VerboseLevel m_vl;  //!< verbose level

  std::string m_mode;

  int m_NinF;             //!< on-site d.o.f.
  int m_Nvol;             //!< volume size
  int m_Ndim;             //!< spacetime dimensions

  AFopr_Domainwall_General<AFIELD>* m_foprdw;

  //! initial setup.
  void init(const Parameters& params);

  //! final tidyup.
  void tidyup();

 public:
  AFopr_OptimalDomainwall(const Parameters& params)
    : AFopr<AFIELD>() { init(params); }

  ~AFopr_OptimalDomainwall() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const double M0, const int Ns,
                const std::vector<int> bc,
                const double b, const double c, 
                const double lambda_min, const double lambda_max);

  void set_config(Field* U) { m_foprdw->set_config(U); }

  void set_config(unique_ptr<Field_G>& U)
  { m_foprdw->set_config(U.get()); }

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(AFIELD& v, const AFIELD& w)
  { m_foprdw->mult(v, w); }

  void mult_dag(AFIELD& v, const AFIELD& w)
  { m_foprdw->mult_dag(v, w); }

  //! mult with specified mode.
  void mult(AFIELD& v, const AFIELD& w, std::string mode)
  { m_foprdw->mult(v, w, mode); }

  void mult_gm5(AFIELD& v, const AFIELD& w)
  { m_foprdw->mult_gm5(v, w); }

  int field_nin(){  return m_NinF; }
  int field_nvol(){ return m_Nvol; }
  int field_nex(){  return m_Ns; }

  //! this returns the number of floating point number operations.
  double flop_count() { return m_foprdw->flop_count(); }

  //! flop-count for specified mode.
  double flop_count(std::string mode)
  { return m_foprdw->flop_count(mode); }

 private:

  void set_optimalDomainwall(const double b, const double c,
                             const double lambda_min,
                             const double lambda_max);

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_OptimalDomainwall(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
                        "OptimalDomainwall", create_object_with_params);
    return init1;
  }
#endif

};
#endif
