/*!
        @file    shiftField_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SHIFTFIELD_EO_INCLUDED
#define SHIFTFIELD_EO_INCLUDED

#include <iostream>

#include "field.h"
#include "index_eo.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift the even-odd field.

/*!
    This class defines the methods that shift a given Field
    instance with even-odd site index in specified direction.
    The field must have half the size of the numbers of sites
    in a node, and field on even(odd) sites are shifted to
    a field on odd (even) sites.
    Thus the method is called with an argument specifying
    one of them: ieo (0: even<-odd, 1: odd<-even).
    The forward shift means, e.g. in mu-direction,
    v(site) = w(site-\hat{mu}), where v is the shifted field
    (output, the first argument) and w the original field
    (input, the second argument).
    At present, there is no implementation for a field which
    has both the even and odd entries (and should be written).
    Names of private functions might be still confusing.
                                   [25 Dec 2011 H.Matsufuru]
   Largely modified toward ver.2.0.
   Multi-threading is incorporated if constructed providing
   Nin.
                                   [27 Nov 2021 H.Matsufuru]
 */

class ShiftField_eo {
 public:
  static const std::string class_name;

 private:
  int m_Nin;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol;
  int m_Nx2, m_Nvol2;

  Index_eo m_index_eo;

  Bridge::VerboseLevel m_vl;

  std::vector<int> m_yzt_eo;

  Field m_wt_x, m_vt_x;  //!< comm. buffer in x-direction
  Field m_wt_y, m_vt_y;  //!< comm. buffer in y-direction
  Field m_wt_z, m_vt_z;  //!< comm. buffer in z-direction
  Field m_wt_t, m_vt_t;  //!< comm. buffer in t-direction

  Field m_we, m_wo;      //!< working field (even/odd)
  Field m_ve, m_vo;      //!< working field (even/odd)
  Field m_w1;            //!< working field (lexical)

 public:

  ShiftField_eo(const int Nin) { init(Nin); }

  ShiftField_eo() { init(); }

 private:
  // non-copyable
  ShiftField_eo(const ShiftField_eo&);
  ShiftField_eo& operator=(const ShiftField_eo&);

 public:
  // forward shift (mu: direction, ieo = 0: even<-odd, 1: odd<-even)
  void forward_h(Field&, const Field&, const int mu, const int ieo);

  // backward shift (mu: direction, ieo = 0: even<-odd, 1: odd<-even)
  void backward_h(Field&, const Field&, const int mu, const int ieo);

  // forward shift with boundary condition
  void forward_h(Field&, const Field&, const int boundary_condition,
                 const int mu, const int ieo);

  // backward shift with boundary condition
  void backward_h(Field&, const Field&, const int boundary_condition,
                  const int mu, const int ieo);

  void forward(Field&, const Field&, const int mu);
  void backward(Field&, const Field&, const int mu);

  void forward(Field&, const Field&,
               const int boundary_condition, const int mu);
  void backward(Field&, const Field&,
                const int boundary_condition, const int mu);

 private:

  void init(const int Nin = 0);

  void up_xh(Field&, const Field&, const int, const int);
  void up_yh(Field&, const Field&, const int, const int);
  void up_zh(Field&, const Field&, const int, const int);
  void up_th(Field&, const Field&, const int, const int);

  void dn_xh(Field&, const Field&, const int, const int);
  void dn_yh(Field&, const Field&, const int, const int);
  void dn_zh(Field&, const Field&, const int, const int);
  void dn_th(Field&, const Field&, const int, const int);
};
#endif
