/*!
      @file    afield_vector-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#include "inline/define_params.h"
#include "inline/define_index.h"

//====================================================================
void afield_set(real_t* RESTRICT v, real_t a, const int nin, const int nst)
{
  int nv = nin * nst;
  for(int iv = 0; iv < nv; ++iv){
    v[iv] = a;
  }

  /*
  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nst; ++ist){
      v[IDXV(nin, nst, in, ist, 0)] = a;
    }
  }
    */

}

//====================================================================
void convert_dev(real_t* RESTRICT v, double* RESTRICT w,
                 int nin, int nvol, int nvol_pad)
{
  int nv1 = nin * nvol;
  int nv2 = nin * nvol_pad;

  for(int ist = 0; ist < nvol_pad; ++ist){
    if(ist < nvol){
      for(int in = 0; in < nin; ++in){
        v[IDXV(nin, nvol_pad, in, ist, 0)] = w[IDX_CORE(nin, in, ist)];
      }
    }else{
      for(int in = 0; in < nin; ++in){
        v[IDXV(nin, nvol_pad, in, ist, 0)] = 0.0;
      }
    }
  }

}

//====================================================================
void reverse_dev(double* RESTRICT v, real_t* RESTRICT w,
                 int nin, int nvol, int nvol_pad)
{
  int nv1 = nin * nvol;
  int nv2 = nin * nvol_pad;

  for(int ist = 0; ist < nvol; ++ist){
    for(int in = 0; in < nin; ++in){
      v[IDX_CORE(nin, in, ist)] = w[IDXV(nin, nvol_pad, in, ist, 0)];
    }
  }

}

//====================================================================
void copy_dev(real_t* RESTRICT v, real_t* RESTRICT w, int nin, int nvol)
{

  int nv = nin * nvol;

  for(int iv = 0; iv < nv; ++iv){
    v[iv] = w[iv];
  }

  /*
  for(int ist = 0; ist < nvol; ++ist){
    for(int in = 0; in < nin; ++in){
      v[IDXV(nin, nvol, in, ist, 0)] = w[IDXV(nin, nvol, in, ist, 0)];
    }
  }
  */
}

//====================================================================
void copy_dev(real_t* RESTRICT v, int nv1, real_t* RESTRICT w, int nv2,
              int nin, int nvol)
{

  int ex1 = nv1/(nin*nvol);
  int ex2 = nv2/(nin*nvol);

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      v[IDXV(nin, nvol, in, ist, ex1)] = w[IDXV(nin, nvol, in, ist, ex2)];
    }
  }

}

//====================================================================
void axpy_dev(real_t *RESTRICT v, int nv1, real_t a,
              real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int nv = nin * nvol;

  int ex1 = nv1/(nin*nvol);
  int ex2 = nv2/(nin*nvol);

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      v[IDXV(nin, nvol, in, ist, ex1)] += a * w[IDXV(nin, nvol, in, ist, ex2)];
    }
  }

}

//====================================================================
void axpy_dev(real_t *RESTRICT v, int nv1, real_t ar, real_t ai,
              real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int nv = nin * nvol;

  int ex1 = nv1/(nin*nvol);
  int ex2 = nv2/(nin*nvol);

  int nin2 = nin/2;

  for(int in = 0; in < nin2; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      real_t wr = w[IDXV(nin, nvol, 2*in,   ist, ex2)];
      real_t wi = w[IDXV(nin, nvol, 2*in+1, ist, ex2)];
      v[IDXV(nin, nvol, 2*in,   ist, ex1)] += ar * wr - ai * wi;
      v[IDXV(nin, nvol, 2*in+1, ist, ex1)] += ai * wr + ar * wi;
    }
  }

}

//====================================================================
void aypx_dev(real_t *RESTRICT v, int nv1, real_t a,
              real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int nv = nin * nvol;

  int ex1 = nv1/(nin*nvol);
  int ex2 = nv2/(nin*nvol);

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      v[IDXV(nin, nvol, in, ist, ex1)]
                              = a * v[IDXV(nin, nvol, in, ist, ex1)]
                                  + w[IDXV(nin, nvol, in, ist, ex2)];
    }
  }

}

//====================================================================
void aypx_dev(real_t *RESTRICT v, int nv1, real_t ar, real_t ai,
              real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int nv = nin * nvol;

  int ex1 = nv1/(nin*nvol);
  int ex2 = nv2/(nin*nvol);

  int nin2 = nin/2;
  for(int in = 0; in < nin2; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      real_t vr = v[IDXV(nin, nvol, 2*in,   ist, ex1)];
      real_t vi = v[IDXV(nin, nvol, 2*in+1, ist, ex1)];
      v[IDXV(nin, nvol, 2*in,   ist,ex1)] = ar * vr - ai * vi
                       + w[nv2 + IDXV(nin, nvol, 2*in,   ist, ex2)];
      v[IDXV(nin, nvol, 2*in+1, ist,ex1)] = ai * vr + ar * vi
                       + w[nv2 + IDXV(nin, nvol, 2*in+1, ist, ex2)];
    }
  }

}

//====================================================================
void scal_dev(real_t* RESTRICT v, int nv1, real_t a, int nin, int nvol)
{
  int nv = nin * nvol;

  int ex1 = nv1/(nin*nvol);

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      v[IDXV(nin, nvol, in, ist, ex1)] *= a;
    }
  }

}

//====================================================================
real_t dot_dev(real_t* RESTRICT v1, real_t* RESTRICT v2,
               int nin, int nvol)
{
  int nv = nin * nvol;

  real_t a = 0.0;

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      a += v1[IDXV(nin, nvol, in, ist, 0)]
                          * v2[IDXV(nin, nvol, in, ist, 0)];
    }
  }

  return a;

}

//====================================================================
real_t norm2_dev(real_t* RESTRICT v1, int nin, int nvol){

  int nv = nin * nvol;

  real_t a = 0.0;

  for(int in = 0; in < nin; ++in){
    for(int ist = 0; ist < nvol; ++ist){
      a += v1[IDXV(nin, nvol, in, ist, 0)] * v1[IDXV(nin, nvol, in, ist, 0)];
    }
  }

  return a;

}

//====================================================================
void dotc_dev(real_t* ar, real_t* ai,
              real_t *RESTRICT v1, real_t *RESTRICT v2,
              int nin, int nst){

  real_t ar2 = 0.0;
  real_t ai2 = 0.0;

  int nv = nin * nst;
  int nin2 = nin/2;

  for(int in = 0; in < nin2; ++in){
    for(int ist = 0; ist < nst; ++ist){
      real_t v1r = v1[IDXV(nin, nst, 2*in,   ist, 0)];
      real_t v1i = v1[IDXV(nin, nst, 2*in+1, ist, 0)];
      real_t v2r = v2[IDXV(nin, nst, 2*in,   ist, 0)];
      real_t v2i = v2[IDXV(nin, nst, 2*in+1, ist, 0)];
      ar2 += v1r * v2r + v1i * v2i;
      ai2 += v1r * v2i - v1i * v2r;
    }

  }

  *ar = ar2;
  *ai = ai2;

}

//====================================================================
void xI_dev(real_t *RESTRICT v1, int nin, int nst){

  int nv = nin * nst;

  int nin2 = nin/2;

  for(int in = 0; in < nin2; ++in){
    for(int ist = 0; ist < nst; ++ist){
      real_t vr = v1[IDXV(nin, nst, 2*in,   ist, 0)];
      real_t vi = v1[IDXV(nin, nst, 2*in+1, ist, 0)];
      v1[IDXV(nin, nst, 2*in,   ist, 0)] = -vi;
      v1[IDXV(nin, nst, 2*in+1, ist, 0)] =  vr;
    }
  }

}

//====================================================================
void conjg_dev(real_t *RESTRICT v1, int nin, int nst){

  int nv = nin * nst;

  int nin2 = nin/2;

  for(int in = 0; in < nin2; ++in){
    for(int ist = 0; ist < nst; ++ist){
      real_t vi = v1[IDXV(nin, nst, 2*in+1, ist, 0)];
      v1[IDXV(nin, nst, 2*in+1, ist, 0)] = -vi;
    }
  }

}

//====================================================================
