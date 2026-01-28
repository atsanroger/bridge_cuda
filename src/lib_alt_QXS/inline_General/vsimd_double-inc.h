/*!
     @file    vsimd_double-inc.h
     @brief
     @author  Hideo Matsufuru (matufuru)
     $LastChangedBy: matufuru $
     @date    $LastChangedDate:: 2023-01-30 10:31:47 #$
     @version $LastChangedRevision: 2465 $
*/

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

typedef unsigned int uint_t;

typedef struct{ double v[VLEND]; } Vsimd_t;

typedef struct{ int v[VLEND]; } Isimd_t;

typedef struct{ unsigned int v[VLEND]; } Usimd_t;

typedef struct{ bool v[VLEND]; } svbool_t;

typedef Vsimd_t svreal_t;

typedef Isimd_t svint_t;

typedef Usimd_t svuint_t;

typedef int int_t;


#endif

