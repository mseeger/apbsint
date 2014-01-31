/* -------------------------------------------------------------------
 * Helper functions for EPTOOLS wrapper functions
 * -------------------------------------------------------------------
 * Macros (this file must be pure C!)
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_HELPER_MACROS_H
#define EPTOOLS_HELPER_MACROS_H

#include <cstdio>
#include <cstring>

// Macros for wrapper functions

// Array argument declarations:
#define W_ARRAY(NAM,TYP) TYP* NAM,int n ## NAM

#define W_DARRAY(NAM) W_ARRAY(NAM,double)

#define W_IARRAY(NAM) W_ARRAY(NAM,int)

#define W_ARR(NAM) NAM,n ## NAM

// Dealing with errors. Requires arguments 'errcode' (int*) and 'errstr'
// (char*), latter a buffer of sufficient size (typically 512)
#define W_ERRCODE errcode

#define W_ERRSTR errstr

#define W_ERRORARGS int* W_ERRCODE,char* W_ERRSTR

#define W_ERRARGS W_ERRCODE,W_ERRSTR

#define W_RETERROR(COD,STR) do {		\
    *W_ERRCODE=COD;				\
    strcpy(W_ERRSTR,STR);			\
    return;					\
  } while (0)

#define W_RETERROR_ARGS(COD,STR,...) do {	\
    *W_ERRCODE=COD;				\
    sprintf(W_ERRSTR,STR,__VA_ARGS__);		\
    return;					\
  } while (0)

#define W_RETOK *W_ERRCODE=0

#define W_CHKSIZE(NAM,SZ,VNAM) do {		\
    if ((n ## NAM) != (SZ))			\
      W_RETERROR(1,VNAM ": Wrong size");	\
  } while (0)

// Flat arrays are often masked by ArrayHandle smart pointers. If the
// array pointer is ?, its length is ?n, the ArrayHandle is ?A
#define W_MASKARRAY(NAM) (NAM ## A).changeRep(NAM,n ## NAM,false)

#endif
