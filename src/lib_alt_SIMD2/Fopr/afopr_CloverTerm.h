/*!
      @file    afopr_CloverTerm.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef AFOPR_CLOVERTERM_H
#define AFOPR_CLOVERTERM_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"

#if defined USE_MPI
#include "lib/Communicator/MPI/communicator_mpi.h"
#include "lib/Communicator/MPI/channel.h"
#else
#include "lib/Communicator/Single/channel.h"
#endif

class Field;
template<typename AFIELD> class AStaple_lex;
template<typename AFIELD> class ASolver;


template<typename AFIELD>
class AFopr_CloverTerm : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd, m_Ndim;
  int m_Nvc, m_Ndf, m_Ndm2;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nstv;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_cSW;               //!< clover coefficient.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;        //!< original gauge config.
  AFIELD m_U;   //!< copied gauge config. with boundary conditions.

  std::string  m_mode;  //!< mult mode

  AStaple_lex<AFIELD> *m_staple;

  ASolver<AFIELD> *m_solver;

  AFIELD m_T, m_Tinv;   //!< clover term (+1) and its inverse.

  AFIELD m_ut1, m_ut2;  //!< working gauge field.

  AFIELD m_v2;          //!< working spinor field.

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  std::vector<Channel*> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

public:
  //! constructor (to be discarded).
  AFopr_CloverTerm(){ init(); }

  //! constructor.
  AFopr_CloverTerm(const Parameters& params){
    init();
    set_parameters(params);
  }

  //! destructor.
  ~AFopr_CloverTerm(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, real_t csW, std::vector<int> bc);

  //! setting gauge configuration.
  void set_config(Field* u);

  void set_config(unique_ptr<Field_G>& u)
  { set_config(u.get()); }

  //! setting fermion boundary condition to gauge field.
  void set_boundary_config(AFIELD& U, const int mu);
    
  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! setting mult mode.
  void set_mode(std::string mode){ m_mode = mode; }

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);
  //  void mult_gm4(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
            const std::string mode);

  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  void mult_csw(AFIELD&, const AFIELD&);
  void multadd_csw(AFIELD&, const AFIELD&);
  void mult_csw_inv(AFIELD&, const AFIELD&);

  //! getting clover term.
  void get_csw(AFIELD& T);

  //! getting the inverse of clover term.
  void get_csw_inv(AFIELD& T);

  //! getting the inverse of clover term with index.
  void get_csw_inv(AFIELD& T, const int);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

 private:
  //! initial setup.
  //  void init(const Parameters& params);
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! setting clover term.
  void set_csw(Field& u);

  //! solve the inverse of clover term.
  void solve_csw_inv();

  //! setting clover term (Dirac repr).
  void set_csw_dirac(Field& u);

  //! setting clover term (Chiral repr).
  void set_csw_chiral(Field& u);

  //! setting field strength.
  void set_fieldstrength(AFIELD& Fst, AFIELD& u, int, int);

  //! multiply diagonal term.
  void multadd_csw(real_t* v2, real_t* v1);

  //! multiply diagonal term at each site (Dirac repr).
  void multadd_csw_dirac(real_t* v2, real_t* v1, int site);

  //! multiply diagonal term at each site (Chiral repr).
  void multadd_csw_chiral(real_t* v2, real_t* v1, int site);

  void mult_gm5(real_t*, real_t*);

};

#endif  // AFOPR_WILSON_H
