/* -------------------------------------------------------------------
 * Helper functions for MEX functions
 * Declarations
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef MEX_HELPER_H
#define MEX_HELPER_H

#include <math.h>
#include "mex.h"
#include "src/eptools/wrap/matrix_types.h"

/*
 * Macros simplifying calling wrapper functions
 * NOTE: M_VARR is for 'fst_vector' type (used in special cases only)
 */

#define M_ARR(NAM) NAM,n ## NAM

#define M_VARR(NAM) NAM.buff,NAM.n

/*
 * M_GET{I|D}{ARRAY|SCAL}: Process input argument
 * - array:        -> int* / double* NAM, size -> n ## NAM
 * - scalar:       -> int/double NAM
 * M_MAKE{I|D}ARRAY: Create return argument (array)
 * - array:        -> int* / double* NAM, size must be in n ## NAM
 *   Buffer created, pointer assigned
 * M_SET{I|D}SCAL: Write return argument (scalar)
 * - scalar:       -> value from NAM
 *
 * NOTE: The int scalar values are double values on the Matlab side, the
 * int32 Matlab type is used for arrays only.
 * NOTE: Every macro starts with ++argidx, so this must be 1 before the
 * intended (e.g. set to -1 for 1st argument)
 */

#define M_GETIARRAY(NAM,STR) n ## NAM = parseInt32Vector(prhs[++argidx],STR,&(NAM),-1)

#define M_GETDARRAY(NAM,STR) n ## NAM = parseDoubleVector(prhs[++argidx],STR,&(NAM),-1)

/*#define M_GETDARRAY(NAM,STR) parseBLASVector(prhs[++argidx],STR,&(NAM),-1,-1)*/

#define M_GETISCAL(NAM,STR) NAM = getScalInt(prhs[++argidx],STR)

#define M_GETDSCAL(NAM,STR) NAM = getScalar(prhs[++argidx],STR)

#define M_MAKEIARRAY(NAM) do {						\
    plhs[++argidx] = mxCreateNumericMatrix(n ## NAM,1,mxINT32_CLASS,mxREAL); \
    NAM = (int*) mxGetData(plhs[argidx]);				\
  } while (0)

#define M_MAKEDARRAY(NAM) do {					\
    plhs[++argidx] = mxCreateDoubleMatrix(n ## NAM,1,mxREAL);	\
    NAM = mxGetPr(plhs[argidx]);				\
  } while (0)

#define M_SETISCAL(NAM) plhs[++argidx] = mxCreateDoubleScalar((double) NAM)

#define M_SETDSCAL(NAM) plhs[++argidx] = mxCreateDoubleScalar(NAM)

/*
 * Exported functions
 *
 * Dealing with last 'errstr' argument: We define macros which set this
 * last argument to 'errstr'. The calling function has to have char* of
 * that name, pointing to the error message buffer.
 */

#define getScalar(arg,name) _getScalar(arg,name,errstr)
double _getScalar(const mxArray* arg,const char* name,char* errstr);

#define getScalInt(arg,name) _getScalInt(arg,name,errstr)
int _getScalInt(const mxArray* arg,const char* name,char* errstr);

#define getVecLen(arg,name) _getVecLen(arg,name,errstr)
int _getVecLen(const mxArray* arg,const char* name,char* errstr);

#define getVecLenAnyType(arg,name) _getVecLenAnyType(arg,name,errstr)
int _getVecLenAnyType(const mxArray* arg,const char* name,char* errstr);

#define getString(arg,name) _getString(arg,name,errstr);
const char* _getString(const mxArray* arg,const char* name,char* errstr);

#define checkMatrix(arg,name,m,n) _checkMatrix(arg,name,m,n,errstr)
void _checkMatrix(const mxArray* arg,const char* name,int m,int n,
		  char* errstr);

#define checkVecIntNonneg(arg,name,n) _checkVecIntNonneg(arg,name,n,errstr)
void _checkVecIntNonneg(const mxArray* arg,const char* name,int n,
			char* errstr);

#define checkVecPosit(arg,name,n) _checkVecPosit(arg,name,n,errstr)
void _checkVecPosit(const mxArray* arg,const char* name,int n,char* errstr);

#define parseBLASMatrix(arg,name,mat,m,n) _parseBLASMatrix(arg,name,mat,m,n,errstr)
void _parseBLASMatrix(const mxArray* arg,const char* name,fst_matrix* mat,
		      int m,int n,char* errstr);

#define parseBLASVector(arg,name,vec,n,bn) _parseBLASVector(arg,name,vec,n,bn,errstr)
void _parseBLASVector(const mxArray* arg,const char* name,fst_vector* vec,
		      int n,int bn,char* errstr);

#define parseDoubleVector(arg,name,vec,n) _parseDoubleVector(arg,name,vec,n,errstr)
int _parseDoubleVector(const mxArray* arg,const char* name,double** vec,int n,
		       char* errstr);

#define parseInt32Vector(arg,name,vec,n) _parseInt32Vector(arg,name,vec,n,errstr)
int _parseInt32Vector(const mxArray* arg,const char* name,int** vec,int n,
		      char* errstr);

/*
 * Allocates void* array of size n and fills with 0. This is used right now
 * to deal with annotation objects, which are not supported in the Matlab
 * interface (for now).
 * NOTE: The array has to be deallocated by 'mxFree' afterwards.
 */
void** getZeroVoidArray(int n);

#endif
