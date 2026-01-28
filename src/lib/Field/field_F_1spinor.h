/*!
        @file    field_F_1spinor.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/


#ifndef FIELD_F_1SPINOR_INCLUDED
#define FIELD_F_1SPINOR_INCLUDED

#include "field_G.h"

#include "Tools/vec_SU_N.h"
using namespace SU_N;
#include "IO/bridgeIO.h"
using Bridge::vout;

//! Staggered-type fermion field.

/*!
  This class defines 1-spinor fermion field, which is mainly used
  by staggered-type fermions.
                                     [28 Dec 2011 H.Matsufuru]

  Added two functions mult_Field_Gn() and mult_Field_Gd() intended to be used
  for fermion gradient flow equation
  \f$\partial_t\chi(t,x)=D_\mu D_\mu\chi(t,x)\f$ since the covariant derivative
  \f$D_\mu\f$ has no spin index.
  [25 May 2017 Y.Taniguchi]
 */
class Field_F_1spinor : public Field
{
 private:
  int m_Nc;   // num of the color elements
  int m_Nc2;  // num of the double color elements

  inline
  int myindex(const int c2, const int site, const int ex) const
  {
    return Field::myindex(c2, site, ex);
  }

 public:

  explicit
  Field_F_1spinor(const int Nvol = CommonParameters::Nvol(), const int Nex = 1) :
    Field(2 * CommonParameters::Nc(), Nvol, Nex),
    m_Nc(CommonParameters::Nc()),
    m_Nc2(2 * CommonParameters::Nc())
  {
  }

  // conversion from Field
  Field_F_1spinor(const Field& x) :
    Field(x),
    m_Nc(CommonParameters::Nc()),
    m_Nc2(2 * CommonParameters::Nc())
  {
  }

  // assignment
  //Field_F_1spinor& operator=(const double a) { field = a; return *this; }
  Field_F_1spinor& operator=(const Field_F_1spinor& v)
  {
    copy(*this, v);
    return *this;
  }

  void reset(int Nvol, int Nex) { Field::reset(m_Nc2, Nvol, Nex); }

  // accessor
  double cmp_r(const int cc, const int site, const int e = 0) const
  {
    return field[myindex(2 * cc, site, e)];
  }

  double cmp_i(const int cc, const int site, const int e = 0) const
  {
    return field[myindex(2 * cc + 1, site, e)];
  }

  void set_r(const int cc, const int site, const int e, const double re)
  {
    field[myindex(2 * cc, site, e)] = re;
  }

  void set_i(const int cc, const int site, const int e, const double im)
  {
    field[myindex(2 * cc + 1, site, e)] = im;
  }

  void set_ri(const int cc, const int site, const int e, const double re, const double im)
  {
    field[myindex(2 * cc, site, e)]     = re;
    field[myindex(2 * cc + 1, site, e)] = im;
  }

  Vec_SU_N vec(const int site, const int e = 0) const
  {
    Vec_SU_N Tmp;

    for (int cc = 0; cc < m_Nc; ++cc) {
      Tmp.set(cc,
              field[myindex(2 * cc, site, e)],
              field[myindex(2 * cc + 1, site, e)]);
    }
    return Tmp;
  }

  void set_vec(const int site, const int e, const Vec_SU_N& F)
  {
    for (int cc = 0; cc < m_Nc; ++cc) {
      field[myindex(2 * cc, site, e)]     = F.r(cc);
      field[myindex(2 * cc + 1, site, e)] = F.i(cc);
    }
  }

  void add_vec(const int site, const int e, const Vec_SU_N& F)
  {
    for (int cc = 0; cc < m_Nc; ++cc) {
      field[myindex(2 * cc, site, e)]     += F.r(cc);
      field[myindex(2 * cc + 1, site, e)] += F.i(cc);
    }
  }

  void clear_vec(const int site, const int e)
  {
    for (int cc = 0; cc < m_Nc2; ++cc) {
      field[myindex(cc, site, e)] = 0.0;
    }
  }
};

// defined as separate functions
//! y=U*x
void mult_Field_Gn(Field_F_1spinor& y, const int ex,
                   const Field_G& u, int ex1,
                   const Field_F_1spinor& x, int ex2);

//! y=U^\dagger*x
void mult_Field_Gd(Field_F_1spinor& y, const int ex,
                   const Field_G& u, int ex1,
                   const Field_F_1spinor& x, int ex2);

#endif
