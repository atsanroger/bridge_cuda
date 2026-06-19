/*!
        @file    shiftField_lex.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SHIFTFIELD_LEX_INCLUDED
#define SHIFTFIELD_LEX_INCLUDED

#include "field.h"
#include "index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift a field in the lexical site index.

/*!
   This class defines the methods which shift a given Field
   instance in the specified direction.
   The forward shift means, e.g. in mu-direction,
   v(site) = w(site-\hat{mu}), where v is the shifted field
   (output, the first argument) and w the original field
   (input, the second argument).
                                   [25 Dec 2011 H.Matsufuru]
   Largely modified toward ver.2.0.
   Multi-threading is incorporated if constructed providing
   Nin.
                                   [27 Nov 2021 H.Matsufuru]
 */

class ShiftField_lex {
 public:
  static const std::string class_name;

 private:
  int m_Nin;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol;

  Bridge::VerboseLevel m_vl;

  Field m_wt_x, m_vt_x;  //!< comm. buffer in x-direction
  Field m_wt_y, m_vt_y;  //!< comm. buffer in y-direction
  Field m_wt_z, m_vt_z;  //!< comm. buffer in z-direction
  Field m_wt_t, m_vt_t;  //!< comm. buffer in t-direction

 public:
  ShiftField_lex() { init(); }

  ShiftField_lex(const int Nin) { init(Nin); }

 private:
  // non-copyable
  ShiftField_lex(const ShiftField_lex&);
  ShiftField_lex& operator=(const ShiftField_lex&);

 public:
  void forward(Field&, const Field&, const int mu);
  void backward(Field&, const Field&, const int mu);

  void forward(Field&, const Field&,
               const int boundary_condition, const int mu);
  void backward(Field&, const Field&,
                const int boundary_condition, const int mu);

 private:

  // initial setup: Nin > 0 leads to resetting the working vectors.
  void init(const int Nin = 0);

  void up_x(Field&, const Field&, const int boundary_condition);
  void up_y(Field&, const Field&, const int boundary_condition);
  void up_z(Field&, const Field&, const int boundary_condition);
  void up_t(Field&, const Field&, const int boundary_condition);

  void dn_x(Field&, const Field&, const int boundary_condition);
  void dn_y(Field&, const Field&, const int boundary_condition);
  void dn_z(Field&, const Field&, const int boundary_condition);
  void dn_t(Field&, const Field&, const int boundary_condition);
};

#endif
