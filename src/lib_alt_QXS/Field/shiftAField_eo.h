/*!
      @file    shiftAField_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef SHIFTAFIELD_EO_INCLUDED
#define SHIFTAFIELD_EO_INCLUDED

#include <vector>

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Methods to shift a field in the even-odd site index.

/*!
   This template class is an alternative to the ShiftField_eo
   class in the core libirary.
   This implementation assumes the QXS site oerdering.
                                   [29 Aug 2021 H.Matsufuru]

   Uses aligned channel.   [5 May 2022 I.anamori]
 */
template<typename AFIELD>
class ShiftAField_eo {
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nin;             //!< internal degree of freedom.
  int m_Nvol, m_Nvol2;
  int m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nx2;
  int m_Nx2v, m_Nyv, m_Nst2v;

  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  std::vector<int> m_Leo;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  ShiftAField_eo(int nin) { init(nin); }

  ShiftAField_eo(int nin, std::vector<int>& bc)
  { init(nin, bc); }

  ~ShiftAField_eo() { tidyup(); }

 private:
  // non-copyable
  ShiftAField_eo(const ShiftAField_eo&);
  ShiftAField_eo& operator=(const ShiftAField_eo&);

 public:

  void forward(AFIELD&, const AFIELD&, const int mu, const int ieo);
  void backward(AFIELD&, const AFIELD&, const int mu, const int ieo);

  void forward(AFIELD&, const int, const AFIELD&, const int,
               const int mu, const int ieo);
  void backward(AFIELD&, const int, const AFIELD&, const int,
                const int mu, const int ieo);

 private:

  //! setup channels for communication.
  void setup_channels();

  void init(int Nin);
  void init(int Nin, std::vector<int>& bc);
  void tidyup();

  void up_xh_simd(real_t *, real_t *, const int);
  void up_yh(real_t *, real_t *, const int);
  void up_zh(real_t *, real_t *, const int);
  void up_th(real_t *, real_t *, const int);
  void dn_xh_simd(real_t *, real_t *, const int);
  void dn_yh(real_t *, real_t *, const int);
  void dn_zh(real_t *, real_t *, const int);
  void dn_th(real_t *, real_t *, const int);

  void up_xh_naive(real_t *, real_t *, const int);
  void up_yh_nv(real_t *, real_t *, const int);
  void up_zh_nv(real_t *, real_t *, const int);
  void up_th_nv(real_t *, real_t *, const int);
  void dn_xh_naive(real_t *, real_t *, const int);
  void dn_yh_nv(real_t *, real_t *, const int);
  void dn_zh_nv(real_t *, real_t *, const int);
  void dn_th_nv(real_t *, real_t *, const int);
};

#endif
