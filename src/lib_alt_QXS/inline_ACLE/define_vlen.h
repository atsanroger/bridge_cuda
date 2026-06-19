/*!
        @file    define_vlen.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/

/*!
    This file provides definitions of VLEN for QXS branch
    of alternative code in Bridge++.
                                 [H.Matsufuru 23 Feb 2021]
*/

#ifndef QXS_DEFINE_VLEN_INCLUDED
#define QXS_DEFINE_VLEN_INCLUDED


#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif

#define VLEND_QXS 8
#define VLENS_QXS 16

#define VLENXD   4
#define VLENYD   2

#define VLENXS   4
#define VLENYS   4

#ifdef QWS_H
// check if consistent with qws.h
#if(VLENS != VLENS_QXS)
#error bad VLENS
#endif
#if(VLEND != VLEND_QXS)
#error bad VLEND
#endif
#undef VLENS
#undef VLEND
#endif

#define VLEND    VLEND_QXS
#define VLENS    VLENS_QXS

#if(VLEND != VLENXD * VLENYD)
#error bad VLENXD * VLENYD
#endif

#if(VLENS != VLENXS * VLENYS)
#error bad VLENXS * VLENYS
#endif

#ifdef USE_QXS_FP16

#ifdef __ARM_FEATURE_SVE
//typedef float16_t half;
typedef _Float16  half;
#else
typedef _Float16  half;
#endif

#define VLENH_QXS 32

#define VLENXH    8
#define VLENYH    4

#define VLENH     VLENH_QXS

#if(VLENH != VLENXH * VLENYH)
#error bad VLENXH * VLENYH
#endif

#endif

#endif
