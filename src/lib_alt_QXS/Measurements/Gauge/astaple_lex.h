/*!
      @file    astaple_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef QXS_ASTAPLE_LEX_INCLUDED
#define QXS_ASTAPLE_LEX_INCLUDED

// core library
#include "lib/Field/field.h"
#include "lib/IO/bridgeIO.h"

// alt-code
#include "lib_alt_QXS/Field/shiftAField_lex.h"

//! Staple construction.

/*!
    This template class is an alternative version of Staple_lex
    class in Bridge++ core library.
                                    [19 Mar 2023 H.Matsufuru]
 */

template<typename AFIELD>
class AStaple_lex
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Ndf, m_Nvol, m_Ndim;

  Bridge::VerboseLevel m_vl;

  ShiftAField_lex<AFIELD> *m_shift;

  AFIELD m_v1, m_v2, m_v3; //!< working fields.
  AFIELD m_Umu;            //!< working fields.

  //! initial setup.
  void init();

  //! final tidy up.
  void tidyup();

 public:

  AStaple_lex(){ init(); }

  ~AStaple_lex(){ tidyup(); }

  //! setting parameters.
  // void set_parameters(const Parameters& params);

  //! calculates plaquette value.
  void plaquette(real_t& plaq, const AFIELD&);

  //! calculates plaquette value.
  void plaquette(real_t& plaq, const Field&);

  //! calculates spatial plaquette value.
  void plaq_s(real_t& plaq_s, const AFIELD&);

  //! calculates temporal plaquette value.
  void plaq_t(real_t& plaq_t, const AFIELD&);

  //! calculates spatial plaquette value.
  void plaq_s_omp(real_t& plaq_s, const AFIELD&);

  //! calculates temporal plaquette value.
  void plaq_t_omp(real_t& plaq_t, const AFIELD&);

  //! calculates spatial plaquette value.
  void plaq_s_impl(real_t& plaq_s, const AFIELD&);

  //! calculates temporal plaquette value.
  void plaq_t_impl(real_t& plaq_t, const AFIELD&);

  //! constructs staple in mu-direction (summing up nu-direction).
  void staple(AFIELD&, const AFIELD&, const int mu);

  //! constructs upper staple in mu-nu plane.
  void upper(AFIELD&, const AFIELD&, const int mu, const int nu);

  //! constructs lower staple in mu-nu plane.
  void lower(AFIELD&, const AFIELD&, const int mu, const int nu);

  //! same as lower, but to be called inside a parallel region.
  void shift_forward(AFIELD& v, const int ex1,
                     const AFIELD& w, const int ex2, const int mu)
  { m_shift->forward(v, ex1, w, ex2, mu); }

  //! same as lower, but to be called inside a parallel region.
  void shift_backward(AFIELD& v, const int ex1,
                     const AFIELD& w, const int ex2, const int mu)
  { m_shift->backward(v, ex1, w, ex2, mu); }

};
#endif
