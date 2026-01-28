/*!
      @file    shiftAField_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef SHIFTAFIELD_LEX_INCLUDED
#define SHIFTAFIELD_LEX_INCLUDED

#include <vector>

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift a field in the lexical site index.

/*!
   This template class is an alternative to the ShiftField
   class in the core libirary.
   This implementation assumes the QXS site oerdering.
                                   [29 Aug 2021 H.Matsufuru]

   Uses aligned channel.   [5 May 2022 I.anamori]
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
  int m_Nxv, m_Nyv, m_Nstv;
  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  ShiftAField_lex(int nin) { init(nin); }

  ShiftAField_lex(int nin, std::vector<int>& bc)
  { init(nin, bc); }

  ~ShiftAField_lex() { tidyup(); }

 private:
  // non-copyable
  ShiftAField_lex(const ShiftAField_lex&);
  ShiftAField_lex& operator=(const ShiftAField_lex&);

 public:
  void forward(AFIELD&, const AFIELD&, const int mu);
  void backward(AFIELD&, const AFIELD&, const int mu);

  void forward(AFIELD&, const int, const AFIELD&, const int,
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

  void up_x_nv(real_t *, real_t *);
  void up_y_nv(real_t *, real_t *);
  void up_z_nv(real_t *, real_t *);
  void up_t_nv(real_t *, real_t *);
  void dn_x_nv(real_t *, real_t *);
  void dn_y_nv(real_t *, real_t *);
  void dn_z_nv(real_t *, real_t *);
  void dn_t_nv(real_t *, real_t *);
};

#endif
