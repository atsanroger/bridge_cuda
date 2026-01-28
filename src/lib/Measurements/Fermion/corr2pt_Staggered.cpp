/*!
        @file    corr2pt_Staggered.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "corr2pt_Staggered.h"

const std::string Corr2pt_Staggered::class_name = "Corr2pt_Staggered";

//====================================================================
void Corr2pt_Staggered::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Corr2pt_Staggered::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Corr2pt_Staggered::meson(std::vector<double>& meson,
                              const std::vector<Field_F_1spinor>& sq1,
                              const std::vector<Field_F_1spinor>& sq2)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = CommonParameters::Nvol();
  const int Lt   = CommonParameters::Lt();
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();

  meson.resize(Lt);

  std::vector<double> corr_local(Nt);
  Field               corrF(2, Nvol, 2);

  {
    const int isrc = 0;
    //for(int isrc = 0; isrc < 2; ++ isrc){

    const int ex = 0;

    for (int site = 0; site < Nvol; ++site) {
      double corr_r = 0.0;
      double corr_i = 0.0;
      for (int c0 = 0; c0 < Nc; ++c0) {
        for (int c1 = 0; c1 < Nc; ++c1) {
          corr_r += sq1[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    + sq1[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_i(c1, site, ex);
          corr_i += sq1[c0 + Nc * isrc].cmp_r(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    - sq1[c0 + Nc * isrc].cmp_i(c1, site, ex)
                    * sq2[c0 + Nc * isrc].cmp_r(c1, site, ex);
        }
      }
      corrF.set(0, site, ex, corr_r);
      corrF.set(1, site, ex, corr_i);
    }

    for (int t = 0; t < Nt; ++t) {
      corr_local[t] = 0.0;

      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          for (int x = 0; x < Nx; ++x) {
            int site = m_index.site(x, y, z, t);
            corr_local[t] += corrF.cmp(0, site, ex);
          }
        }
      }
    }


    const int ipe_t = Communicator::ipe(3);

    std::vector<double> corr_global(Lt);
    for (int t = 0; t < Lt; ++t) {
      corr_global[t] = 0.0;
    }
    for (int t = 0; t < Nt; ++t) {
      int t2 = t + ipe_t * Nt;
      corr_global[t2] = corr_local[t];
    }
    for (int t = 0; t < Lt; ++t) {
      double cr_r = corr_global[t];

      meson[t] = Communicator::reduce_sum(cr_r);
    }
  }
}


//====================================================================
//============================================================END=====
