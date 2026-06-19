/*!
        @file    index_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef INDEX_EO_INCLUDED
#define INDEX_EO_INCLUDED

#include <vector>

#include "Field/index_lex.h"
#include "Field/field.h"

#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Even-odd site index.

/*!
    This class defines even-odd site index.
    Only the site degree of freedom is concerned.
    Nx (x-extent inside a node) must be even in the present
    implementation.
    Coverting from and reverting to the lexical site index are
    implemented as member functions of this class.
    In present implementation, there is no superclass structure,
    and thus polymorphism is not available.
    Some of method names might be confusing; restructuring may
    be helpful.
                                      [25 Dec 2011 H.Matsufuru]
    Multi-threading is applied to methods.
                                      [29 Nov 2021 H.Matsufuru]
*/
class Index_eo {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;
  int m_Nx2, m_Nvol2;
  int m_node_eo;   //!< {0, 1}: local origin is on even/odd side

  std::vector<int> m_yzt_eo;
  std::vector<int> m_site_up;
  std::vector<int> m_site_dn;

  Index_lex m_index_lex;

 public:

  Index_eo() { init(); }

  int leo(const int y, const int z, const int t) const
  { return m_yzt_eo[y + m_Ny * (z + m_Nz * t)]; }

  int site(const int x2, const int y, const int z, const int t,
           const int ieo) const
  { return x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo; }

  int site(const int is, const int ieo) const
  { return is + m_Nvol2 * ieo; }

  int site_up(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_up[s] + m_Nvol2 * (1 - ieo);
  }

  int site_xup(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_up[s] + m_Nvol2 * (1 - ieo);
  }

  int site_dn(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_dn[s] + m_Nvol2 * (1 - ieo);
  }

  int site_xdn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_dn[s] + m_Nvol2 * (1 - ieo);
  }

  int siteh(const int x2, const int y, const int z, const int t)
  const
  { return x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)); }

  int siteh_up(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_up[s];
  }

  int siteh_xup(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_up[s];
  }

  int siteh_dn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_dn[s];
  }

  int siteh_xdn(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
    return m_site_dn[s];
  }

  //! initial setup.
  void init();

  void convertField(Field& eo, const Field& lex);
  void convertField(Field& eo, const Field& lex, const int ieo);

  void reverseField(Field& lex, const Field& eo);
  void reverseField(Field& lex, const Field& eo, const int ieo);

  void splitField(Field& e, Field& o, const Field& eo);

  void mergeField(Field& eo, const Field& e, const Field& o);
};
#endif
