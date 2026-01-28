/*!
      @file    afield_dd_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2024-03-21 21:09:15 #$
      @version $LastChangedRevision: 2587 $
*/


// Note that
//  * ieo = 0 or 1 -> operation only on even or odd blocks
//  * ieo < 0 (typically -1) -> operation on all the blocks
// The parity of origin block must be handled in the parent function.

//====================================================================
void block_axpy_eo(real_t *restrict vp, real_t *restrict out,
                   real_t *restrict wp, real_t fac,
                   int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  // In this function, out is an array of size (2 * Nblock),
  // assuming Nblock complex values.

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  int nv = Nin * Nst * Nex;

#pragma acc data present(vp[0:nv], out[0:2*Nblock], wp[0:nv])	\
                 copyin(Nx, Ny, Nz, Nt, Nst, ieo, fac,	\
                        Bx, By, Bz, Bt, Bsize, NBx, NBy, NBz, Nblock)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
    for(int site2 = 0; site2 < Nst; ++site2){
      int bsite = site2 % Bsize;
      int block = site2/Bsize;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);

      int ieo_skip = 1 - ieo;
      int jeo = (ibx + iby + ibz + ibt) % 2;
      if(jeo == ieo_skip) continue;

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite/Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt/By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt/Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      real_t ar = out[2*block]   * fac;
      real_t ai = out[2*block+1] * fac;

      for(int ex = 0; ex < Nex; ++ex){
        int nv1 = Nin * Nst * ex;
        int nv2 = Nin * Nst * ex;
        for(int in2 = 0; in2 < Nin/2; ++in2){
          int inr = 2 * in2;
          int ini = 2 * in2 + 1;
          real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
          real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
          real_t wtr = wp[nv2 + IDX2(Nin, inr, site)];
          real_t wti = wp[nv2 + IDX2(Nin, ini, site)];
          vtr += ar * wtr - ai * wti;
          vti += ar * wti + ai * wtr;
          vp[nv1 + IDX2(Nin, inr, site)] = vtr;
          vp[nv1 + IDX2(Nin, ini, site)] = vti;

        }
      }

    }

  } // acc parallel

}

//====================================================================
void block_scal_eo(real_t *restrict vp, real_t *restrict out,
                   int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  // In this function, out is an array of size (2 * Nblock).

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  int nv = Nin * Nst * Nex;

#pragma acc data present(vp[0:nv], out[0:2*Nblock]) \
                 copyin(Nx, Ny, Nz, Nt, Nst, ieo, \
                        Bx, By, Bz, Bt, Bsize, NBx, NBy, NBz, Nblock)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
    for(int site2 = 0; site2 < Nst; ++site2){
      int bsite = site2 % Bsize;
      int block = site2/Bsize;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);

      int ieo_skip = 1 - ieo;
      int jeo = (ibx + iby + ibz + ibt) % 2;
      if(jeo == ieo_skip) continue;

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite/Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt/By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt/Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      real_t ar = out[2*block];
      real_t ai = out[2*block+1];

      for(int ex = 0; ex < Nex; ++ex){
        int nv1 = Nin * Nst * ex;
        for(int in2 = 0; in2 < Nin/2; ++in2){
          int inr = 2 * in2;
          int ini = 2 * in2 + 1;
          real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
          real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
          vp[nv1 + IDX2(Nin, inr, site)] = ar * vtr - ai * vti;
          vp[nv1 + IDX2(Nin, ini, site)] = ar * vti + ai * vtr;
        }
      }

    }

  } // acc parallel

}

//====================================================================
void block_dotc_eo(real_t *restrict out,
                   real_t *restrict vp, real_t *restrict wp,
                   int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  // In this function, out is an array of size (2 * Nst).

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  int nv = Nin * Nst * Nex;

#pragma acc data present(out[0:2*Nst], vp[0:nv], wp[0:nv])	\
                 copyin(Nx, Ny, Nz, Nt, Nst, ieo, \
                        Bx, By, Bz, Bt, Bsize, NBx, NBy, NBz, Nblock)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
    for(int site2 = 0; site2 < Nst; ++site2){
      int bsite = site2 % Bsize;
      int block = site2/Bsize;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);

      int ieo_skip = 1 - ieo;
      int jeo = (ibx + iby + ibz + ibt) % 2;
      if(jeo == ieo_skip) continue;

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite/Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt/By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt/Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      real_t xtr = 0.0;
      real_t xti = 0.0;

      for(int ex = 0; ex < Nex; ++ex){
        int nv1 = Nin * Nst * ex;
        int nv2 = Nin * Nst * ex;
        for(int in2 = 0; in2 < Nin/2; ++in2){
          int inr = 2 * in2;
          int ini = 2 * in2 + 1;
          real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
          real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
          real_t wtr = wp[nv2 + IDX2(Nin, inr, site)];
          real_t wti = wp[nv2 + IDX2(Nin, ini, site)];
          xtr += vtr * wtr + vti * wti;
          xti += vtr * wti - vti * wtr;
        }
      }

      out[2 * site2]     = xtr;
      out[2 * site2 + 1] = xti;

    }

  } // acc parallel

}

//====================================================================
void block_norm2_eo(real_t *restrict out, real_t *restrict vp, 
                   int Nin, int Nex,
                   int *Nsize, int *block_size, int ieo)
{
  // In this function, out is an array of size (2 * Nst).

  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

  int nv = Nin * Nst * Nex;

#pragma acc data present(out[0:2*Nst], vp[0:nv])  \
                 copyin(Nx, Ny, Nz, Nt, Nst, ieo, \
                        Bx, By, Bz, Bt, Bsize, NBx, NBy, NBz, Nblock)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
    for(int site2 = 0; site2 < Nst; ++site2){
      int bsite = site2 % Bsize;
      int block = site2/Bsize;

      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);

      int ieo_skip = 1 - ieo;
      int jeo = (ibx + iby + ibz + ibt) % 2;
      if(jeo == ieo_skip) continue;

      int kx   = bsite % Bx;
      int ix   = kx + Bx * ibx;
      int kyzt = bsite/Bx;
      int ky   = kyzt % By;
      int iy   = ky + By * iby;
      int kzt  = kyzt/By;
      int kz   = kzt % Bz;
      int iz   = kz + Bz * ibz;
      int kt   = kzt/Bz;
      int it   = kt + Bt * ibt;
      int site = ix + Nx * (iy + Ny * (iz + Nz * it));

      real_t xtr = 0.0;

      for(int ex = 0; ex < Nex; ++ex){
        int nv1 = Nin * Nst * ex;
        for(int in2 = 0; in2 < Nin/2; ++in2){
          int inr = 2 * in2;
          int ini = 2 * in2 + 1;
          real_t vtr = vp[nv1 + IDX2(Nin, inr, site)];
          real_t vti = vp[nv1 + IDX2(Nin, ini, site)];
          xtr += vtr * vtr + vti * vti;
        }
      }

      out[2 * site2]     = xtr;
      out[2 * site2 + 1] = 0.0;

    }

  } // acc parallel

}

//====================================================================
void reduce_block(real_t *restrict out, real_t *restrict out_red,
		  // int Nblock, int Nst)
                  int* Nsize, int* block_size, int ieo)
{
  int Nx = Nsize[0];
  int Ny = Nsize[1];
  int Nz = Nsize[2];
  int Nt = Nsize[3];
  int Nst  = Nx * Ny * Nz * Nt;

  // size of block
  int Bx = block_size[0];
  int By = block_size[1];
  int Bz = block_size[2];
  int Bt = block_size[3];
  int Bsize = Bx * By * Bz * Bt;

  // numbers of blocks
  int NBx  = Nsize[0]/block_size[0];
  int NBy  = Nsize[1]/block_size[1];
  int NBz  = Nsize[2]/block_size[2];
  int NBt  = Nsize[3]/block_size[3];
  int Nblock = NBx * NBy * NBz * NBt;

#pragma acc data present(out[0:2*Nblock], out_red[0:2*Nst]) \
                 copyin(Nblock, Bsize, NBx, NBy, NBz, ieo)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {

#pragma acc loop
    for(int block = 0; block < Nblock; ++block){
      int ibx = block % NBx;
      int iby = (block/NBx) % NBy;
      int ibz = (block/(NBx * NBy)) % NBz;
      int ibt = block/(NBx * NBy * NBz);

      int ieo_skip = 1 - ieo;
      int jeo = (ibx + iby + ibz + ibt) % 2;

      if(jeo == ieo_skip){
        out[2*block]   = 0.0;
        out[2*block+1] = 0.0;
      }else{
        real_t vtr = 0.0;
        real_t vti = 0.0;
        for(int bsite = 0; bsite < Bsize; ++ bsite){
          real_t *vp = &out_red[2 * Bsize * block];
          vtr += vp[2 * bsite];
          vti += vp[2 * bsite +1];
        }
        out[2*block]   = vtr;
        out[2*block+1] = vti;
      }
    }
  }

}

//============================================================END=====
