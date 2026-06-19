/*!
      @file    afield_copy_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-11-05 14:28:33 #$
      @version $LastChangedRevision: 2618 $
*/

//====================================================================
void copy(double* restrict v, int nv1, float* restrict w, int nv2,
              int nin, int nvol)
{

  int nv = nin * nvol;

#pragma acc data present(v[nv1:nv], w[nv2:nv]) \
                 copyin(nin, nvol, nv1, nv2)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void copy(float* restrict v, int nv1, double* restrict w, int nv2,
              int nin, int nvol)
{

  int nv = nin * nvol;

#pragma acc data present(v[nv1:nv], w[nv2:nv]) \
                 copyin(nin, nvol, nv1, nv2)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
