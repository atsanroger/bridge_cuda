/*!
        @file    vsimd_float-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
	$LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-01-30 10:31:47 #$
        @version $LastChangedRevision: 2465 $
*/

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

typedef unsigned int uint_t;

typedef struct{ float v[VLENS]; } Vsimd_t;

typedef struct{ int v[VLENS]; } Isimd_t;

typedef struct{ unsigned int v[VLENS]; } Usimd_t;

typedef struct{ bool v[VLENS]; } svbool_t;

typedef Vsimd_t svreal_t;

typedef Isimd_t svint_t;

typedef Usimd_t svuint_t;

typedef int int_t;


#endif
