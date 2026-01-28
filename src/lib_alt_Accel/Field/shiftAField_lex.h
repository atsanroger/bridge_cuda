/*!
        @file    shiftAField_lex.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#ifndef ACCEL_SHIFTAFIELD_LEX_INCLUDED
#define ACCEL_SHIFTAFIELD_LEX_INCLUDED

#include <vector>

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift a field in the lexical site index.

/*!
   This template class is an alternative version of ShiftField
   in the Bridge core libirary.
   This implementation assumes the complex array which is
   integrated in the inner-most SIMD vector.
   Differently from the Field version, the Nin size must be
   specified at the construction.
   Note that the Nin/2 is placed as the size of second degree
   of freedom if the field array.
   The boundary condition is also given at the construction.
                                   [25 Dec 2017 H.Matsufuru]
 */
template<typename AFIELD>
class ShiftAField_lex {

 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nin;             //!< internal degree of freedom.
  int m_Nvol;
  int m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nxv, m_Nstv;
  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];
  int m_bc[4];
  int m_bc2[4];

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  ShiftAField_lex(int nin) { init(nin); }

  ShiftAField_lex(int nin, std::vector<int>& bc)
  { init(nin, bc); }

  ~ShiftAField_lex()
  { tidyup(); }

 private:
  // non-copyable
  ShiftAField_lex(const ShiftAField_lex&);
  ShiftAField_lex& operator=(const ShiftAField_lex&);

 public:
  void forward(AFIELD&, const AFIELD&, const int mu);
  void backward(AFIELD&, const AFIELD&, const int mu);

  void forward( AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);
  void backward(AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);

 private:

  //! setup channels for communication.
  void setup_channels();

  void init(int Nin);
  void init(int Nin, std::vector<int>& bc);

  void tidyup();

  void up_x(real_t *, real_t *);
  void up_y(real_t *, real_t *);
  void up_z(real_t *, real_t *);
  void up_t(real_t *, real_t *);
  void dn_x(real_t *, real_t *);
  void dn_y(real_t *, real_t *);
  void dn_z(real_t *, real_t *);
  void dn_t(real_t *, real_t *);

};

#endif
