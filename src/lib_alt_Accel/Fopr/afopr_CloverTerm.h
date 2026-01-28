/*!
      @file    afopr_CloverTerm.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-09-10 17:39:07 #$
      @version $LastChangedRevision: 2537 $
*/

#ifndef ACCEL_AFOPR_CLOVERTERM_H
#define ACCEL_AFOPR_CLOVERTERM_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

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
  enum Gamma_repr {DIRAC, CHIRAL} m_repr;  //!< gamma matrix representation

  int m_Nc, m_Nd, m_Ndim;
  int m_Nvc, m_Ndf, m_Ndm2;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_cSW;               //!< clover coefficient.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  //  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;        //!< original gauge config.
  AFIELD m_U;   //!< copied gauge config. with boundary conditions.

  std::string  m_mode;  //!< mult mode

  AStaple_lex<AFIELD> *m_staple;
  ASolver<AFIELD> *m_solver;

  AFIELD m_T, m_Tinv;   //!< clover term (+1) and its inverse.

  AFIELD m_v1, m_v2;    //!< working spinor fields.

  AFIELD m_ut1, m_ut2;  //!< working gauge fields.

  AFIELD m_U2, m_F2;    //!< working gauge fields.

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];

  //  std::vector<int> m_Nbdsize;
  //  using Channel = Channel_impl<std::allocator<real_t> >;
  //  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
//  ChannelSet chset_send, chset_recv;

public:
  //! constructor (to be discarded).
  //  AFopr_CloverTerm(){ init(); }

  //! constructor.
  AFopr_CloverTerm(const Parameters& params){
    init(params);
  }

  //! destructor.
  ~AFopr_CloverTerm(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, real_t csW, std::vector<int> bc);

  //! setting gauge configuration.
  void set_config(Field* u);

  //  void set_config(unique_ptr<Field_G>& u)
  //  { set_config(u.get()); }

  //! setting fermion boundary condition to gauge field.
  void set_boundary_config(AFIELD& U, const int mu);
    
  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! setting mult mode.
  void set_mode(std::string mode);

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
  void init(const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  void set_config_omp(Field* u);
  void set_config_impl(Field* u);

  //! setting clover term.
  void set_csw(Field& u);

  //! solve the inverse of clover term.
  void solve_csw_inv();
  void solve_csw_inv_dirac();
  void solve_csw_inv_chiral();

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


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(
                                      const Parameters& params)
  { return new AFopr_CloverTerm(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register(
                         "CloverTerm", create_object_with_params);
    return init1;
  }
#endif
  
};

#endif  // ACCEL_AFOPR_CLOVERTERM_H
