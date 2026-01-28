/*!
      @file    prefetch.h
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

/*
 Copyright
   Bridge++ project and RIKEN (2022)

 Licence: GPL
    see README.txt and LICENSE for more details
*/
#pragma once

#if !defined(DISABLE_PREFETCH)

#define __prefetch_load_luinv(a, offset)                                   \
  {                                                                        \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 0)], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 4)], 0, 2); \
  }
#define __prefetch_write_luinv(a, offset)                                  \
  {                                                                        \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 0)], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 4)], 1, 2); \
  }

#define __prefetch_load_luinv_l1(a, offset)                                \
  {                                                                        \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 0)], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 4)], 0, 3); \
  }
#define __prefetch_write_luinv_l1(a, offset)                               \
  {                                                                        \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 0)], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (site + 1) + VLEN * (offset + 4)], 1, 3); \
  }

#define __prefetch_load_hop_u_l2(a, dir, idx)                                       \
  {                                                                                 \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 0 * 64], 0, 2); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 1 * 64], 0, 2); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 2 * 64], 0, 2); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 3 * 64], 0, 2); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 4 * 64], 0, 2); \
  }
#define __prefetch_load_hop_vec_l2(a, idx, is)                       \
  {                                                                  \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 0 * 64], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 1 * 64], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 2 * 64], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 3 * 64], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 4 * 64], 0, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 5 * 64], 0, 2); \
  }

#define __prefetch_write_hop_vec_l2(a, idx, is)                      \
  {                                                                  \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 0 * 64], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 1 * 64], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 2 * 64], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 3 * 64], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 4 * 64], 1, 2); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 5 * 64], 1, 2); \
  }

#define __prefetch_load_hop_u_l1(a, dir, idx)                                       \
  {                                                                                 \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 0 * 64], 0, 3); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 1 * 64], 0, 3); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 2 * 64], 0, 3); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 3 * 64], 0, 3); \
    __builtin_prefetch(&a[NDF * Nst2 * (dir) + VLEN * NDF * (idx) + 4 * 64], 0, 3); \
  }
#define __prefetch_load_hop_vec_l1(a, idx, is)                       \
  {                                                                  \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 0 * 64], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 1 * 64], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 2 * 64], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 3 * 64], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 4 * 64], 0, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 5 * 64], 0, 3); \
  }

#define __prefetch_write_hop_vec_l1(a, idx, is)                      \
  {                                                                  \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 0 * 64], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 1 * 64], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 2 * 64], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 3 * 64], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 4 * 64], 1, 3); \
    __builtin_prefetch(&a[Nin5 * (idx) + Nin4 * is + 5 * 64], 1, 3); \
  }

#define __prefetch_load_hop2_buf_x_l2(a, idx, is, skip)     \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 2); \
  }

#define __prefetch_load_hop2_buf_y_l2(a, idx, is, skip)     \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 2); \
  }

#define __prefetch_load_hop2_buf_zt_l2(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 2); \
    __builtin_prefetch(&a[idx + skip * is + 1 * 64], 0, 2); \
    __builtin_prefetch(&a[idx + skip * is + 2 * 64], 0, 2); \
  }

#define __prefetch_load_hop2_buf_x_l1(a, idx, is, skip)     \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 3); \
  }

#define __prefetch_load_hop2_buf_y_l1(a, idx, is, skip)     \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 3); \
  }

#define __prefetch_load_hop2_buf_zt_l1(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 0, 3); \
    __builtin_prefetch(&a[idx + skip * is + 1 * 64], 0, 3); \
    __builtin_prefetch(&a[idx + skip * is + 2 * 64], 0, 3); \
  }

#define __prefetch_write_hop1_buf_x_l2(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 2); \
  }

#define __prefetch_write_hop1_buf_y_l2(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 2); \
  }

#define __prefetch_write_hop1_buf_zt_l2(a, idx, is, skip)   \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 2); \
    __builtin_prefetch(&a[idx + skip * is + 1 * 64], 1, 2); \
    __builtin_prefetch(&a[idx + skip * is + 2 * 64], 1, 2); \
  }

#define __prefetch_write_hop1_buf_x_l1(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 3); \
  }

#define __prefetch_write_hop1_buf_y_l1(a, idx, is, skip)    \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 3); \
  }

#define __prefetch_write_hop1_buf_zt_l1(a, idx, is, skip)   \
  {                                                         \
    __builtin_prefetch(&a[idx + skip * is + 0 * 64], 1, 3); \
    __builtin_prefetch(&a[idx + skip * is + 1 * 64], 1, 3); \
    __builtin_prefetch(&a[idx + skip * is + 2 * 64], 1, 3); \
  }

#else

#define __prefetch_load_luinv(a, offset)
#define __prefetch_write_luinv(a, offset)
#define __prefetch_load_luinv_l1(a, offset)
#define __prefetch_write_luinv_l1(a, offset)
#define __prefetch_load_hop_u_l2(a, dir, idx)
#define __prefetch_load_hop_vec_l2(a, idx, is)
#define __prefetch_write_hop_vec_l2(a, idx, is)
#define __prefetch_load_hop_u_l1(a, dir, idx)
#define __prefetch_load_hop_vec_l1(a, idx, is)
#define __prefetch_write_hop_vec_l1(a, idx, is)
#define __prefetch_load_hop2_buf_x_l2(a, skip, idx, is)
#define __prefetch_load_hop2_buf_y_l2(a, idx, is, skip)
#define __prefetch_load_hop2_buf_zt_l2(a, idx, is, skip)
#define __prefetch_load_hop2_buf_x_l1(a, idx, is, skip)
#define __prefetch_load_hop2_buf_y_l1(a, idx, is, skip)
#define __prefetch_load_hop2_buf_zt_l1(a, idx, is, skip)
#define __prefetch_write_hop1_buf_x_l2(a, idx, is, skip)
#define __prefetch_write_hop1_buf_y_l2(a, idx, is, skip)
#define __prefetch_write_hop1_buf_zt_l2(a, idx, is, skip)
#define __prefetch_write_hop1_buf_x_l1(a, idx, is, skip)
#define __prefetch_write_hop1_buf_y_l1(a, idx, is, skip)
#define __prefetch_write_hop1_buf_zt_l1(a, idx, is, skip)
#endif
