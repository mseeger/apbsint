//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class NumberFormats
 * ------------------------------------------------------------------- */

#ifndef NUMBERFORMATS_H
#define NUMBERFORMATS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header
#include "lhotse/FileUtils.h"

/**
 * Template class allowing type-checked access to 'FileUtils', to store/
 * load sequences of elementary types in a machine-independent way.
 * <p>
 * NOTE: All (binary) file I/O of LHOTSE should use exclusively these
 * methods, except in very special situations in which large temporary
 * copies have to be dumped which are guaranteed to be read on the same
 * architecture and NOT be stored for later use!
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
template<class T> class NumberFormats
{
public:
  // Public static methods

  /**
   * Stores sequence of length 'n', starting from 'data', using a
   * machine-indep. format.
   * If 'index' is given, then 'n' is the length of the index, and 'data' is
   * indexed in the sense that element i is found at data[index[i]]. If 'index'
   * is not given, then 'step' is the step size for the vector buffer, i.e.
   * element i is at data+i*step. 'step' can be negative, but not ==0.
   * <p>
   * NOTE: The default implementation here just calls 'saveSeq', passing
   * sizeof(T). Override this by specialisation if this is not appropriate for
   * T.
   *
   * @param os    Output file stream
   * @param data  Pointer to sequence
   * @param n     Length of sequence. Def.: 1
   * @param step  Step size. Ignored if 'index' is given. Def.: 1
   * @param index See above. Def.: 0
   * @param fsize Byte size of T in file. Def.: -1 (system byte size)
   */
  static void save(ofstream& os,const T* data,int n=1,int step=1,
		   const int* index=0,int fsize=-1) {
    FileUtils::saveSeq(os,(const char*) data,sizeof(T),n,step,index,fsize);
  }

  /**
   * Loads sequence of length 'n', writing to 'data', using a machine-indep.
   * format.
   * If 'index' is given, then 'n' is the length of the index, and 'data' is
   * indexed in the sense that element i is found at data[index[i]]. If 'index'
   * is not given, then 'step' is the step size for the vector buffer, i.e.
   * element i is at data+i*step. 'step' can be negative, but not ==0.
   * <p>
   * NOTE: The default implementation here just calls 'loadSeq', passing
   * sizeof(T). Override this by specialisation if this is not appropriate for
   * T.
   *
   * @param is    Input file stream
   * @param data  Pointer to sequence
   * @param n     Length of sequence. Def.: 1
   * @param step  Step size. Ignored if 'index' is given. Def.: 1
   * @param index See above. Def.: 0
   * @param fsize Byte size of T in file. Def.: -1 (system byte size)
   */
  static void load(ifstream& is,T* data,int n=1,int step=1,const int* index=0,
		   int fsize=-1)
  {
    FileUtils::loadSeq(is,(char*) data,sizeof(T),n,step,index,fsize);
  }
};

#endif
