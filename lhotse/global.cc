//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Global definitions
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#ifdef MATLAB_MEX
#include "lhotse/matif/mex_for_cpp.h"
#endif

// See AC_F77_DUMMY_MAIN in the autoconf manual
#ifdef F77_DUMMY_MAIN
extern "C" int F77_DUMMY_MAIN() { return 1; }
#endif

// Messages to stdout

void printMsgStdout(const char* msg)
{
#ifdef MATLAB_MEX
  //PRINTGLOBALLOG(msg);
  mexPrintf(msg); mexPrintf("\n");
#else
  cout << msg << endl;
#endif
}

#if defined(USE_OWN_MEMMAN) || (defined(MATLAB_MEX) && defined(USE_MATLAB_MM))

// Global overload of new/delete and their array variants.

void* operator new(std::size_t n) throw(std::bad_alloc)
{
#ifdef USE_OWN_MEMMAN
  if (FixedMemManager::isMMInit())
    return FixedMemManager::alloc(n);
#endif
#if !defined(MATLAB_MEX) || !defined(USE_MATLAB_MM)
  return malloc(n);
#else
  void* ptr=mxMalloc(n);
  mexMakeMemoryPersistent(ptr);
  return ptr;
#endif
}

void* operator new[](std::size_t n) throw(std::bad_alloc)
{
#ifdef USE_OWN_MEMMAN
  if (FixedMemManager::isMMInit())
    return FixedMemManager::alloc(n);
#endif
#if !defined(MATLAB_MEX) || !defined(USE_MATLAB_MM)
  return malloc(n);
#else
  void* ptr=mxMalloc(n);
  mexMakeMemoryPersistent(ptr);
  return ptr;
#endif
}

void operator delete(void* ptr) throw()
{
#ifdef USE_OWN_MEMMAN
  if (FixedMemManager::isMMInit())
    FixedMemManager::dealloc(ptr);
#endif
#if !defined(MATLAB_MEX) || !defined(USE_MATLAB_MM)
  free(ptr);
#else
  mxFree(ptr);
#endif
}

void operator delete[](void* ptr) throw()
{
#ifdef USE_OWN_MEMMAN
  if (FixedMemManager::isMMInit())
    FixedMemManager::dealloc(ptr);
#endif
#if !defined(MATLAB_MEX) || !defined(USE_MATLAB_MM)
  free(ptr);
#else
  mxFree(ptr);
#endif
}

#endif

// DEBUG code to find memory leaks

#ifdef DEBUG_TRACKHANDLES

MAP_TYPE(int,debugType) debugMem;

void debugMemPrintStats()
{
  int len=0,i,j,sz=14;
  int num[sz];
  long total[sz];
  bool tagged;

  for (i=0; i<sz; i++) {
    num[i]=0; total[i]=0;
  }
  printMsgStdout("*** debugMemPrintStats:");
  char msg[100];
  MAP_CONSTITER(int,debugType) it=debugMem.begin();
  while (it!=debugMem.end()) {
    len++;
    i=(int) (*it).second.tag;
    if (i>=sz)
      throw InternalException(EXCEPT_MSG(""));
    num[i]++;
    total[i]+=(long) (*it).second.sz;
    // Tag==8 is of interest:
    if (i==8) {
      const int* cnt=(*it).second.cnt;
      MemWatchBase* ptr=(*it).second.ptr;
      sprintf(msg,"    Tag 8: %d %d %d %d %d %d %d [%d,sz=%d]",cnt[0],cnt[1],
	      cnt[2],cnt[3],cnt[4],cnt[5],cnt[6],ptr->getRefCount(),
	      (*it).second.sz);
      printMsgStdout(msg);
    }
    ++it;
  }
  printMsgStdout("*** Summary:");
  for (i=0; i<sz; i++) {
    j=(num[i]>0)?((int) (total[i]/((long) num[i]))):0;
    sprintf(msg,"    %d: %d (avg. sz. %d)",i,num[i],j);
    printMsgStdout(msg);
  }
}

#endif
