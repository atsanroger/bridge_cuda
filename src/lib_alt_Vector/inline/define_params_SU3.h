/*!
        @file    define_params_SU3.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#ifndef DEFINE_PARAMS_SU3_INCLUDED
#define DEFINE_PARAMS_SU3_INCLUDED

// lattice parameters assuming gauge group of SU(3)
#define  NC      3
#define  NVC     6
#define  NDF    18
#define  NDF2    9
#define  ND      4
#define  ND2     2
#define  NCD    12
#define  NVCD   24
#define  NDIM    4

// vector length for SX-Aurora
#define  VLEN  256
#define  VPAD    16
//#define  VSHIFT   128
#define  VSHIFT   0

// OpenACC parameters: to be dicarded
//#define  NWP     1
#define  NWP     16

#define NUM_WORKERS 1
#define VECTOR_LENGTH 32

#endif

//============================================================END=====
