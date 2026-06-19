/*!
        @file    bridge_init_factory_alt_SIMD2.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$

        @version $LastChangedRevision: 2160 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_SIMD2_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_SIMD2_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_SIMD2();

#endif /* USE_FACTORY_AUTOREGISTER */
  
#ifdef DEBUG
void bridge_report_factory_alt_SIMD2();
#endif

#endif /* USE_FACTORY */


#endif
