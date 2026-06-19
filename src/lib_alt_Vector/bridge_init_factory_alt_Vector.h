/*!
        @file    bridge_init_factory_alt_Vector.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$

        @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_VECTOR_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_VECTOR_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_Vector();

#endif /* USE_FACTORY_AUTOREGISTER */
  
#ifdef DEBUG
void bridge_report_factory_alt_Vector();
#endif

#endif /* USE_FACTORY */


#endif
