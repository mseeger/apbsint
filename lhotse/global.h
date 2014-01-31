//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Global header
 * ------------------------------------------------------------------- */

#ifndef GLOBAL_H
#define GLOBAL_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

/**
 * Provides declarations of global classes and other global declarations.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */

/*
 * Things related to config.h macros provided by 'configure'. This is done
 * only if HAVE_CONFIG is defined.
 * LHOTSE should also compile stand-alone, without 'configure' support.
 * To this end, if HAVE_CONFIG is not given, required config.h macros
 * not of HAVE_$$$ form are given default values here if they are not
 * supplied. In general, use the -D compiler option to supply macros.
 *
 * The following config.h macros are presently used:
 * - HAVE_CONFIG: In every source file, config.h is included if this is
 *   defined
 * - HAVE_DEBUG: LHOTSE is run in debug mode. This has a number of
 *   consequences. HAVE_DEBUG should not be set when compiling a stable
 *   application
 * - HAVE_LIBGSL: The GNU scientific library is available.
 *   Note that this is somewhat obsolete, since the GSL is required in
 *   general
 * - HAVE_INLINE: Compiler does inlining. This flag is used by the GSL
 *   headers. This may depend on the optimization level of the compiler.
 *   Use -DHAVE_INLINE compiler option to define this one
 * - HAVE_OS_???: Flag indicating the operating system for which LHOTSE
 *   is compiled. See below
 * - HAVE_FORTRAN: If this is not defined, external Fortran code is not
 *   used, corr. wrappers throw a 'NotImplemException'.
 *   NOTE: This is about Fortran source code shipped with LHOTSE, not
 *   about library functions from BLAS, LAPACK, without which LHOTSE does
 *   not compile at all
 * - F77_FUNC: A Fortran function 'foo' is called by
 *     F77_FUNC(foo,FOO) (...)
 */

/*
 * Principal flags
 *
 * Flags marked with (*) are serviced through the Makefile. They should
 * not be defined in general if the code is compiled using this Makefile.
 *
 * - MATLAB_MEX (*):
 *   Flag must be set iff current program is to become a Matlab MEX script.
 *   Calls to MEX interface functions should be wrapped.
 *   NOTE: Define USE_MATLAB_MM as well if the standard C++ memory manager
 *   is to be replaced by Matlab's MM.
 *   NOTE: If MATLAB_MEX is set, ther DYNCAST macro translates to a static
 *   rather than a dynamic cast, because dynamic_cast does not work for
 *   MEX files. Code relying on dynamic_cast for type checking will not
 *   work, it has to be wrapped and dealt with specially.
 * - USE_MATLAB_MM (*):
 *   Only considered if MATLAB_MEX is defined. If set, the new/delete
 *   operators are redefined to use the Matlab memory manager, with all
 *   allocations being made persistent. Also, a non-standard allocator is
 *   used for STL objects if these are created using the macros XXX_TYPE
 *   macros defined below.
 *   NOTE: Following the Matlab documentation, one should always use the
 *   Matlab MM within MEX programs. But then, memory allocated with the
 *   standard C++ manager seems to remain persistent as well, and using the
 *   Matlab MM for everything seems to be instable.
 *   If the MEX routine does not maintain persistent memory, the flag should
 *   always be set.
 * - MATLAB_VER65:
 *   Only considered if MATLAB_MEX is defined. Set this variable iff you are
 *   using Matlab 6.5 or later. Set by default.
 * - USE_OWN_MEMMAN:
 *   If this is defined, we replace the global memory manager by our own
 *   which is more efficient if a large number of allocations/deallocations
 *   come at a small number of fixed small sizes (see 'FixedMemManager').
 *   ATTENTION: Does not work properly right now!
 * - MATLAB_DEBUG / MATLAB_DEBUG_OLD:
 *   MatlabDebug code supposed to be active should be enclosed in
 *   MATLAB_DEBUG. Old MatlabDebug code not to be used should be enclosed in
 *   MATLAB_DEBUG_OLD.
 *   MATLAB_DEBUG_OLD must NEVER be defined! MATLAB_DEBUG should be defined
 *   by default, and undefined only if no MatlabDebug code should be active
 *   at all.
 * - MATLABDEBUG_USEMEX
 *   If this and MATLAB_MEX are set, the 'MatlabDebug' services are run by
 *   calling back Matlab directly using MEX commands s.a. 'mexEvalString'.
 *   This is much more efficient than the file-based mode, but less robust
 *   and only applicable to MEX functions (normal code can be wrapped into
 *   a dummy MEX function to use this mode)
 * - CHECKASSERTS:
 *   Assertion checking, see 'AssertMethod'.
 *   Assertions are checked only if CHECKASSERTS is defined. In this case,
 *   CHECKASSERTS is a status flag:
 *   - 0: MYASS just prints message if assertion fails
 *   - 1: MYASS prints message and terminates program [default]
 *   ATTENTION: At present, none of the LHOTSE classes are child of
 *   'AssertMethod' or implement 'assertm'.
 *   NOTE: MYASS should be used in time-critical, internal code. CHECKASSERTS
 *   can be removed after debugging, to get rid of these checks.
 *   ==> Right now, CHECKASSERTS is defined (==1) iff HAVE_DEBUG is defined.
 * - NAMESPACE: Set if namespace facility available. If NAMESPACE is not
 *   set, name clashes are not allowed.
 *   NOTE: This is from the early design and has never been used, may not be
 *   consequently implemented! It requires further code modifications in order
 *   to support namespaces. LHOTSE at present does not use namespaces.
 * - DEBUG_TRACKHANDLES:
 *   Debug code to track memory regions used by ArrayHandle and the BaseXXX
 *   matrix/vector classes. Lots of stuff!
 */

#ifdef HAVE_DEBUG
#define CHECKASSERTS 1
#endif
//#define MATLAB_MEX
//#define USE_MATLAB_MM
#define MATLAB_VER65
//#define USE_OWN_MEMMAN // DO NOT USE RIGHT NOW!
#define MATLAB_DEBUG
#define MATLABDEBUG_USEMEX
//#define NAMESPACE
//#define DEBUG_TRACKHANDLES

/*
 * Compiling on different systems:
 * LHOTSE is supported under Linux only right now. It is tested and does
 * properly compile only for this case (HAVE_OS_LINUX).
 * Still, things which are system-dependent and we are aware of, can make
 * use of the constants here. The default definition here is ignored if
 * any of them is set already (by the Makefile).
 */
#if !defined(HAVE_OS_LINUX) && !defined(HAVE_OS_WINDOWS) && !defined(HAVE_OS_SUN_SOLARIS) && !defined(HAVE_OS_ALPHA_OSF)
#  define HAVE_OS_LINUX
#endif

// See AC_F77_DUMMY_MAIN in the autoconf manual
// Definition is in global.cc
#ifdef F77_DUMMY_MAIN
extern "C" int F77_DUMMY_MAIN();
#endif

// If F77_FUNC is not defined in config.h, we provide a variant which
// is OK for Linux or Windows
#ifndef F77_FUNC
#ifdef HAVE_OS_LINUX
#  define F77_FUNC(name,cname) name ## _
#else
#  define F77_FUNC(name,cname) name
#endif
#endif

// Frequently used modules from STL

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <climits>
#include <cfloat>

/*
 * The std namespace
 * Starting from GCC 3, the "std" namespace is taken literally. We could just
 * write
 *   using namespace std;
 * but this is not recommended. We deal with this by the following "using"
 * clauses for the most frequently used names. All others have to be
 * prefixed explicitly with "std::".
 */

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::istream;
using std::pair;
using std::ofstream;
using std::ostream;
using std::string;
using std::stringstream;
using std::less;
using std::greater;
using std::equal_to;
using std::not_equal_to;
using std::less_equal;
using std::greater_equal;
using std::logical_and;
using std::logical_or;
using std::logical_not;
using std::plus;
using std::minus;
using std::multiplies;
using std::divides;
using std::modulus;
using std::negate;
using std::bind1st;
using std::bind2nd;
using std::ptr_fun;
using std::not1;
using std::not2;

/*
 * Global namespace macros. If NAMESPACE is defined, namespaces are used.
 * Namespaces cannot be handled by older compilers, and newer compilers
 * are often buggy w.r.t. namespaces. Naming clashes should be avoided
 * until namespaces are fully established! However, code should be developed
 * using namespaces, because it can be integrated into a larger project
 * much easier then.
 * - use "NS(nsname)ident..." instead of "nsname::ident..."
 * - use "USING(nsname);" instead of "using namespace nsname;"
 * - use "BEGINNS(nsname)" instead of "namespace nsname {"
 * - use "ENDNS" instead of the "}" corr. to a "BEGINNS(nsname)"
 * Note: Emacs automatic identing is spoiled by BEGINNS/ENDNS. It is wise
 * to comment these statements during editing and uncomment before
 * compilation. If NAMESPACE is not set, you can leave the comments in
 */

#ifdef NAMESPACE
#define NS(nsname) nsname::
#define USING(nsname) using namespace nsname
#define BEGINNS(nsname) namespace nsname {
#define ENDNS }
#else
#define NS(nsname)
#define USING(nsname)
#define BEGINNS(nsname)
#define ENDNS
#endif

/*
 * STL-related types:
 * If code is used in a MEX script, all STL objects which should remain
 * persistent beyond the return of the MEX function have to use the Matlab MM
 * and make alloc. memory persistent. This is done by providing the
 * non-standard allocator 'MatlabAlloc'. Use types/macros defined below, the
 * correct allocator is selected based on MATLAB_MEX.
 *
 * IMPORTANT: Use these types and macros exclusively in code which may be used
 * in MEX files. Extend macro/typedef list whenever new STL classes are
 * required. The original STL classes can be used for purely local objects
 * which are not class members.
 * - types 'my_XXX': for non-template STL classes (e.g. string)
 * - macros: for template STL classes
 *
 * ATTENTION: Use of "typename" keyword:
 * If you get a lot of "implicit typename is deprecated" warnings during
 * compilation, and possibly the compiler stops due to lack of memory, here
 * is the reason:
 * Suppose you write within a template definition (with param. T) something
 * like:
 *   VEC_CONSTITER(T) it=...;
 * The true type is obtained from a member of std::vector<T>. To simplify
 * the compiler work, the C++ standard requires you to explicitly tell that
 * to the compiler:
 *   typename VEC_CONSTITER(T) it=...;
 * ==> "typename XXX" means that XXX has to be considered as a typename.
 * Confusingly, "typename" must NOT be used if you are not in a template
 * context (don't ask me why!), so
 *   typename VEC_CONSTITER(int) it=...;
 * gives you a "using `typename' outside of template" error.
 * NOTE: This affects the XXX_CONSTITER and XXX_ITER macros, but NOT the
 * XXX_TYPE macros for which "typename" must not be used.
 */
#if defined(MATLAB_MEX) && defined(USE_MATLAB_MM)

template<class T> class MatlabAlloc;
typedef std::basic_string<char,std::char_traits<char>,MatlabAlloc<char> > my_string;

#define MAP_CONSTITER(key,data) std::map<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >::const_iterator
#define MAP_ITER(key,data) std::map<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >::iterator
#define MAP_TYPE(key,data) std::map<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >
#define MMAP_CONSTITER(key,data) std::multimap<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >::const_iterator
#define MMAP_ITER(key,data) std::multimap<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >::iterator
#define MMAP_TYPE(key,data) std::multimap<key,data,std::less<key >,MatlabAlloc<pair<const key,data > > >
#define VEC_CONSTITER(data) std::vector<data,MatlabAlloc<data > >::const_iterator
#define VEC_ITER(data) std::vector<data,MatlabAlloc<data > >::iterator
#define VEC_TYPE(data) std::vector<data,MatlabAlloc<data > >
#define SET_CONSTITER(data) std::set<data,std::less<data >,MatlabAlloc<data > >::const_iterator
#define SET_ITER(data) std::set<data,std::less<data >,MatlabAlloc<data > >::iterator
#define SET_TYPE(data) std::set<data,std::less<data >,MatlabAlloc<data > >
#define LIST_CONSTITER(data) std::list<data,MatlabAlloc<data > >::const_iterator
#define LIST_ITER(data) std::list<data,MatlabAlloc<data > >::iterator
#define LIST_TYPE(data) std::list<data,MatlabAlloc<data > >

#else

typedef std::string my_string;

#define MAP_CONSTITER(key,data) std::map<key,data,std::less<key > >::const_iterator
#define MAP_ITER(key,data) std::map<key,data,std::less<key > >::iterator
#define MAP_TYPE(key,data) std::map<key,data,std::less<key > >
#define MMAP_CONSTITER(key,data) std::multimap<key,data,std::less<key > >::const_iterator
#define MMAP_ITER(key,data) std::multimap<key,data,std::less<key > >::iterator
#define MMAP_TYPE(key,data) std::multimap<key,data,std::less<key > >
#define VEC_CONSTITER(data) std::vector<data >::const_iterator
#define VEC_ITER(data) std::vector<data >::iterator
#define VEC_TYPE(data) std::vector<data >
#define SET_CONSTITER(data) std::set<data >::const_iterator
#define SET_ITER(data) std::set<data >::iterator
#define SET_TYPE(data) std::set<data >
#define LIST_CONSTITER(data) std::list<data >::const_iterator
#define LIST_ITER(data) std::list<data >::iterator
#define LIST_TYPE(data) std::list<data >

#endif

// Global macros

/*
 * Dynamic cast: Due to a bug, MEX scripts crash when dynamic_cast is used.
 * We use static cast if MATLAB_MEX is set.
 * ==> ATTENTION: Does not solve the problem! Using static instead of
 *     dynamic cast is often NOT appropriate and never safe!
 * ==> If A virtual base of B, the compiler rejects static cast from
 *     A* to B*
 */
#ifdef MATLAB_MEX
#define DYNCAST(type,arg) ((type*) arg)
#else
#define DYNCAST(type,arg) (dynamic_cast<type*>(arg))
#endif

/*
 * Throw exceptions with file name and line number. Usage:
 *   throw ExceptionName(EXCEPT_MSG("Message"));
 */
#define EXCEPT_MSG(msg) msg,__FILE__,__LINE__

/*
 * 'close' on a 'ifstream' / 'ofstream' does not reset the stream
 * state. To properly close a stream which can then be reused, 'clear'
 * has to be called as well (one of the most reported "bugs" for GCC 3!).
 * Use CLOSESTR.
 */
#define CLOSESTR(st) do { (st).close(); (st).clear(); } while(0)

// Global typedefs

typedef unsigned char uchar;
typedef unsigned int  uint;
typedef unsigned long ulong;

/*
 * Global functions
 * - printMsgStdout: Messages to stdout (the usual 'cout << ...' is not
 *   appropriate with MATLAB_MEX)
 */

extern void printMsgStdout(const char*);

#ifdef DEBUG_TRACKHANDLES

/*
 * Tracking down memory leaks:
 * There is some code available to do this:
 * - DEBUG_TRACKHANDLES: See below
 * - FIXEDMEMMANAGER_COLLECT_STATS in 'FixedMemManager': Collects stats
 *   on the usage of the MM for fixed-sized objects. Can detect build-up
 *   of small memory blocks which are not dealloc.
 */

/* DEBUG code to track handling of memory regions via 'MemWatcher'
   (see 'ArrayHandle' code). Useful to find memory leaks in this system.
   Global map 'debugMem': Ptr. mem. region to debug structure (debugType).
   Memory region receives tag which identifies its creator. Tag==0 means
   the creator is unknown. Tag values can be passed to constructor of
   'MemWatcher'. Current setting (CHANGES!):
   0:  Unknown, or 'ArrayHandle'
   1:  BaseVector::resizeSave
   2:  BaseVector::expand
   3:  BaseVector::insert(1)
   4:  BaseVector::insert(2)
   5:  BaseVector::init
   6:  BaseVector::ensureCapacity
   7:  BaseMatrix constructor(rows,cols)
   8:  BaseMatrix::ensureCapacity
   9:  BaseMatrix::resizeSave
   10: BaseMatrix::expand
   11: BaseMatrix::shrinkBuffer
   12: BaseMatrix constructor(rows,cols,a)
   13: BaseMatrix::resize
   A second mechanism tracks causes of the ref. counter (MemWatchBase)
   being incr./decr. This is only active for regions with specific tag
   value. First, identify tag value for problematic regions. Then, switch
   on cause tracking for this tag value.
   In the moment (CHANGES!): 7.
   Uses 'debugCause' member in 'MatTimeStamp'. Possible causes are:
   0: Unknown
   1: ArrayHandle methods
   2: Matrix masking part of matrix
   3: Vector masking column of matrix
   4: Vector masking row of matrix
   5: Vector masking diagonal of matrix
   6: Vector masking part of vector
   Each entry in 'debugMem' keeps a counter for each of the causes. Counters
   are increased together with pcount in watcher, decreased when the cause
   object deletes its association with the buffer. An object which curr. does
   not have any others masking parts of it, should have all counters ==0.
   NOTE: Cause counters are updated only for specific tag values (see above).
*/

class MemWatchBase;

struct debugType {
  int sz;
  uchar tag;
  MemWatchBase* ptr;
  int cnt[7];

  debugType() {
    sz=0; tag=0; ptr=0;
    for (int i=0; i<7; i++) cnt[i]=0;
  }

  debugType(const debugType& a) {
    for (int i=0; i<7; i++) cnt[i]=a.cnt[i];
    tag=a.tag; sz=a.sz; ptr=a.ptr;
  }
};
extern MAP_TYPE(int,debugType) debugMem;
extern void debugMemPrintStats();

#endif

/*
 * Memory manager 'FixedMemManager' for fixed-sized objects:
 * The MM is maintained in the static members/methods of 'FixedMemManager'.
 * It initializes itself upon first usage. We globally overload new/delete
 * (see global.cc, global_mem.h) to use this MM once properly initialized.
 * NOTE: Done iff USE_OWN_MEMMAN is defined.
 * NOTE: Does not work properly in the moment. Do not use!
 */

// LHOTSE global include files

#include "lhotse/StandardException.h" // Superclass of LHOTSE exceptions
#include "lhotse/exceptions.h"        // Global exceptions
#if defined(MATLAB_MEX) && defined(USE_MATLAB_MM)
#include "lhotse/matif/MatlabAlloc.h" // Matlab memory allocator
#endif
#ifdef USE_OWN_MEMMAN
#include "lhotse/FixedMemManager.h"   // MM for fixed-sized objects
#endif
#include "lhotse/global_mem.h"        // Global overload of new/delete (code
                                      // in global.cc)
#include "lhotse/AssertMethod.h"      // Assertion method, checking
#include "lhotse/Handle.h"            // Handles (smart pointers)
#include "lhotse/ArrayHandle.h"       // Handles for arrays, mem. watchers
#include "lhotse/ArrayPtrHandle.h"    // Handles for inhomogenous arrays
#include "lhotse/LogFile.h"           // Logfile support
#include "lhotse/DefaultLogs.h"       // Default logfiles (global,error)
#ifdef HAVE_DEBUG
#include "lhotse/DebugVars.h"         // Debug facilities
#endif
#include "lhotse/Range.h"             // Ranges
#include "lhotse/Interval.h"          // Intervals
#include "lhotse/FuncObjects.h"       // Global function objects, predicates
#include "lhotse/NullaryFunc.h"       // Nullary functions
#include "lhotse/AccumulFunc.h"       // Accumulators

#endif
