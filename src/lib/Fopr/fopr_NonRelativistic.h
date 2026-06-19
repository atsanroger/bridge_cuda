/*!
      @File    fopr_NonRelativistic.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef FOPR_NONRELATIVISTIC_INCLUDED
#define FOPR_NONRELATIVISTIC_INCLUDED

#include "Fopr/fopr.h"
#include "Field/field_G.h"

//! Nonrelativistic QCD fermion operator.

/*!
  This class implements the nonrelativistic QCD (NRQCD)
  fermion operator.
  References of implementation:
   [1] K.I. Ishikawa et al., Phys. Rev. D 56 (1997) 7028.
   [2] S. Hashimoto et al.,  Phys. Rev. D 58 (1998) 014502.
  Convention of correction terms is based on [1].
  - evolution type-A corresponds to [1]
  - evolution type-B corresponds to [2]
  We thank Takafumi Aoki for checking the code.
                         [02 Jul 2021/30 Dec 2022 H.Matsufuru]
 */

class Fopr_NonRelativistic : public Fopr
{
 public:
  static const std::string class_name;

 private:
  double m_MQ;                    //!< heavy quark mass
  int m_nstab;                    //!< stabilization parameter
  int m_num_correct;              //!< number of correction terms
  double m_u0;                    //!< mean-field value

  Bridge::VerboseLevel m_vl;      //!< verbose level

  std::vector<double> m_coeff;    //!< coefficients of correction terms
  std::vector<int> m_boundary;    //!< boundary conditions
  std::string m_evolution_type;   //!< type of evolution equation
  std::string m_correction_terms; //!< order of correction terms

  std::string m_repr;             //!<  gamma matrix repr. (set to Dirac)

  std::string m_mode;             //!<  mode of multiplication

  int m_Nc, m_Nvc, m_Nd, m_Nd2;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nvol, m_Ndim;
  int m_Nspc;                         //!< spatial volume (m_Nx * m_Ny * m_Nz)
  int m_Ltime;                        //!< global temporal extent

  int m_itime_src;                    //!< source time slice

  Field_G m_U;                        //!< gauge configuration

  std::vector<Field_G> m_Fstr;        //!< Field strenth (0-2: B, 3-5: E)

  std::vector<double> buf1_x, buf2_x; //!< communication buffer in x-dir.
  std::vector<double> buf1_y, buf2_y; //!< communication buffer in y-dir.
  std::vector<double> buf1_z, buf2_z; //!< communication buffer in z-dir.

  std::vector<Field> m_Xt;            //!< heavy quark field.
  Field m_Wt1, m_Wt2, m_Wt3;          //!< working field on each time slice.
  Field m_Yt1, m_Yt2, m_Yt3;          //!< working field on each time slice.
  Field m_Zt1, m_Zt2;                 //!< working field on each time slice.

 public:

  //! constructor without argument

  /*
  Fopr_NonRelativistic(): Fopr()
  { init(); }
  */

  //! constructor with parameter
  Fopr_NonRelativistic(const Parameters& params) : Fopr()
  { init(params); }

  //! destructor
  ~Fopr_NonRelativistic() { tidyup(); }

  //! seting parameters with a Parameter object
  void set_parameters(const Parameters& params);

  //! seting parameters with individual values
  void set_parameters(const double MQ,
                      const int nstab,
                      const double m_u0,
                      const std::vector<double> coeff,
                      const std::vector<int> bc,
                      const std::string evolution_type,
                      const std::string correction_terms);

  void get_parameters(Parameters& params) const;

  //! setting gauge configuration and field strength
  void set_config(Field *U);

  void set_config(unique_ptr<Field_G>& U) { set_config(U.get()); }

  //! setting mult mode: 'Evolve' and 'Rotation'
  void set_mode(const std::string mode);

  std::string get_mode() const { return m_mode; }

  //! mult with m_mode
  void mult(Field& v, const Field& w);

  //! not available
  void mult_dag(Field& v, const Field& w);

  //! mult with give mode
  void mult(Field& v, const Field& w, const std::string mode);

  //! not available
  void mult_dag(Field& v, const Field& w, const std::string mode);

  //! multiply $\gamma_5$ : not available
  void mult_gm5(Field& v, const Field& w);

  //! transporter in upper direction : not available
  void mult_up(const int mu, Field& v, const Field& w);

  //! transporter in lower direction : not available
  void mult_dn(const int mu, Field& v, const Field& w);

  //! volume size as a 4D operator
  int field_nvol() { return m_Nvol; }

  //! inner degree of freedom as a 4-spinor
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! extra degree of freedom
  int field_nex() { return 1; }

  //! number of floating point operations: not implemented
  double flop_count();

  //! flop_count for given mode: not implemented
  double flop_count(const std::string mode);

 private:

  //! initial setup.
  void init(const Parameters& params);

  //! final clean up.
  void tidyup();

  //! evolution equation (facade).
  void evolve(Field&, const Field&);

  //! extract source time slice and 3D source vector from 4D source.
  void set_source(Field& Xt, const Field& b);

  //!  called from evolve() and switchs evlution equation
  void evolve_impl(Field&, const Field&, const int itime);

  //! evolution equation according to [1]
  void evolve_typeA(Field&, const Field&, const int itime);

  //! evolution equation according to [2]
  void evolve_typeB(Field&, const Field&, const int itime);

  //! evolution with kinetic term
  void evolve_H0(Field&, const Field&, const int itime);

  //! evolution with one time slice ahead
  void evolve_U4(Field&, const Field&, const int itime);

  //! evolution with correction terms
  void evolve_deltaH(Field&, const Field&, const int itime);

  //! correction term deltaH(1): spin-magnetic interaction
  void add_deltaH1(Field&, const Field&, const int itime);

  //! correction term deltaH(2)
  void add_deltaH2(Field&, const Field&, const int itime);

  //! correction term deltaH(3)
  void add_deltaH3(Field&, const Field&, const int itime);

  //! correction term deltaH(5)
  void add_deltaH5(Field&, const Field&, const int itime);

  //! correction terms deltaH(4) + deltaH(6)
  void add_deltaH46(Field&, const Field&, const int itime);

  //! second order covariant derivative
  void calc_Delta2(Field&, const Field&, const int itime);

  //! first order covariant derivative
  void calc_Delta1(Field&, const Field&, const int idir, const int itime);

  //! fourth order covariant derivative
  void calc_Delta4(Field&, const Field&, const int itime);

  //! field rotation (jd = 1: normal, -1: conjugate)
  void rotation(Field& v, const Field& w, const int jd);

  void add_R2(Field& v, const Field& w, const int jd, const int itime);

  void add_R3(Field& v, const Field& w, const int itime);

  void add_R4(Field& v, const Field& w, const int itime);

  void add_R5(Field& v, const Field& w, const int itime);

  //! multiplying Pauli matrix $\sigma_1$
  void mult_sigma1(Field& Xt, const Field& Yt);

  //! multiplying Pauli matrix $\sigma_2$
  void mult_sigma2(Field& Xt, const Field& Yt);

  //! multiplying Pauli matrix $\sigma_3$
  void mult_sigma3(Field& Xt, const Field& Yt);

  //! shifting 3D field forward in idir-direction
  void shift_forward(Field& Xt, const Field& Yt,
                     const int idir, const int itime);

  //! shifting 3D field backward in idir-direction
  void shift_backward(Field& Xt, const Field& Yt,
                      const int idir, const int itime);

  //! implemetation of shifting field
  void shift_xup(Field& Xt, const Field& Yt, const int itime);
  void shift_xdn(Field& Xt, const Field& Yt, const int itime);
  void shift_yup(Field& Xt, const Field& Yt, const int itime);
  void shift_ydn(Field& Xt, const Field& Yt, const int itime);
  void shift_zup(Field& Xt, const Field& Yt, const int itime);
  void shift_zdn(Field& Xt, const Field& Yt, const int itime);

  //! nultiplication of gauge field
  void mult_Gn(Field& Xt, const Field& Yt, const int idir, const int itime);

  //! nultiplication of gauge field (hermitian conjugate)
  void mult_Gd(Field& Xt, const Field& Yt, const int idir, const int itime);

  //! nultiplication of field strength
  void mult_F(Field& Xt, const Field& Yt, const int icomp, const int itime);


#ifdef USE_FACTORY
 private:
  static Fopr *create_object_with_params(const Parameters& params)
  { return new Fopr_NonRelativistic(params); }

 public:
  static bool register_factory()
  {
    bool init = Fopr::Factory_params::Register("NonRelativistic",
                                               create_object_with_params);
    return init;
  }
#endif
};
#endif
