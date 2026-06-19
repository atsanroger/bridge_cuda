/*!
        @file    afopr_Chebyshev.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-07-05 20:48:53 #$
        @version $LastChangedRevision: 2530 $
*/

#ifndef AFOPR_CHEBYSHEV_INCLUDED
#define AFOPR_CHEBYSHEV_INCLUDED

#include "lib/Fopr/afopr.h"

#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

#include <vector>
#include <string>
using std::string;

class Field;

//! Chebyshev polynomial.

/*!
    This class implements the chebyshev polynomial of a fermion
    operator. Based on Fopr_Chebyshev and modified so that
    double/float can be swtched.
                                     [12 Jul 2018 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Chebyshev : public AFopr<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  int m_NinF, m_Nvol, m_NexF;
  int m_Npcb;
  real_t m_Vthrs, m_Vmax;
  real_t m_Fcb1, m_Fcb2;

  AFopr<AFIELD> *m_fopr;

  std::string m_mode;         //!< mult mode

  Bridge::VerboseLevel m_vl;  //!< verbose level

  std::vector<AFIELD> m_dj;

 public:
  //! standard constructor.
  AFopr_Chebyshev(AFopr<AFIELD> *fopr,
                  const Parameters& params)
    : AFopr<AFIELD>(), m_fopr(fopr)
  { init(params); }

  DEPRECATED
  AFopr_Chebyshev(AFopr<AFIELD> *fopr)
    : AFopr<AFIELD>(), m_fopr(fopr)
  { init(); }

  ~AFopr_Chebyshev() { tidyup(); }

  //! setting parameters with a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters with values.
  void set_parameters(const int Np,
                      const real_t v_threshold,
                      const real_t v_max);

  void get_parameters(Parameters& params) const;

  //! setting gauge configuration.
  void set_config(Field *U);

  void set_mode(std::string mode);

  std::string get_mode() const;

  void mult(AFIELD& v, const AFIELD& w);

  void mult_dag(AFIELD& v, const AFIELD& w)
  {
    mult(v, w);
  }

  void mult(real_t& v, const real_t x);

  real_t mult(const real_t x);

  //! returns true if additional field conversion is needed.
  virtual bool needs_convert()
  { return m_fopr->needs_convert(); }

  //! converts a Field object into other format if necessary.
  virtual void convert(AFIELD& v, const Field& w)
  { m_fopr->convert(v, w); }

  //! reverses to a Field object from other format if necessary.
  virtual void reverse(Field& v, const AFIELD& w)
  { m_fopr->reverse(v, w); }

  int field_nin()  { return m_NinF; }
  int field_nvol() { return m_Nvol; }
  int field_nex()  { return m_NexF; }

  double flop_count();

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! initial setup (obsolete).
  void init();

  //! final cleanup.
  void tidyup();

#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object(AFopr<AFIELD> *fopr)
  { return new AFopr_Chebyshev<AFIELD>(fopr); }

  static AFopr<AFIELD> *create_object_with_params(AFopr<AFIELD> *fopr,
                                                  const Parameters& params)
  { return new AFopr_Chebyshev<AFIELD>(fopr, params); }

 public:
  static bool register_factory()
  {
    bool                init = true;
    init &= AFopr<AFIELD>::Factory_fopr::Register(
      "Chevyshev", create_object);
    init &= AFopr<AFIELD>::Factory_fopr_params::Register(
      "Chevyshev", create_object_with_params);
    return init;
  }
#endif
};

#endif
