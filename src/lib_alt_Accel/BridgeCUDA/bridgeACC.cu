/*!
      @file    bridgeACC.cu
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-07-23 18:49:38 #$
      @version $LastChangedRevision: 2534 $
*/


#include <stdio.h>
#include <algorithm>
#include "bridgeACC.h"

namespace BridgeACC {

void*  ptr_host[MAX_DEVPTR_REF];
void*  ptr_dev[MAX_DEVPTR_REF];
size_t ptr_size[MAX_DEVPTR_REF];
int devptr_max = 0;

//====================================================================
void* dev_ptr(void* phost)
{
/*
  for(int i = 0; i < devptr_max; ++i){
    if(ptr_host[i] == phost) return ptr_dev[i];
  }
  printf("no registered device pointer: %p\n", phost);
  exit(1);
  return 0;
*/

  for(int i = 0; i < devptr_max; ++i){
    char* upper = ((char*)ptr_host[i] + ptr_size[i]);
    char* lower = (char*)ptr_host[i];
    if((char*)phost >= lower && (char*)phost < upper){
      size_t offset = (char*)phost - (char*)ptr_host[i];
      return (void*)((char*)ptr_dev[i] + offset);
    }
  }
  printf("no registered device pointer: %p\n", phost);
  exit(1);
  return 0;


}

}
//============================================================END=====
