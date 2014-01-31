//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Declarations for new/delete overload
 * ------------------------------------------------------------------- */

#ifndef GLOBALMEM_H
#define GLOBALMEM_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#if defined(USE_OWN_MEMMAN) || (defined(MATLAB_MEX) && defined(USE_MATLAB_MM))

extern void* operator new(std::size_t n) throw(std::bad_alloc);
extern void operator delete(void* ptr) throw();
extern void* operator new[](std::size_t n) throw(std::bad_alloc);
extern void operator delete[](void* ptr) throw();

#endif

#endif
