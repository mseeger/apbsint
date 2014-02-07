/* -------------------------------------------------------------------
 * Helper functions for MEX functions
 * Implementation
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "matlab/mex/mex_helper.h"

/*
 * Exported functions
 */

double _getScalar(const mxArray* arg,const char* name,char* errstr)
{
  if (!mxIsDouble(arg) || mxGetM(arg)!=1 || mxGetN(arg)!=1) {
    sprintf(errstr,"Expect double scalar for %s",name);
    mexErrMsgTxt(errstr);
  }
  return *mxGetPr(arg);
}

int _getScalInt(const mxArray* arg,const char* name,char* errstr)
{
  double val,temp;

  if (!mxIsDouble(arg) || mxGetM(arg)!=1 || mxGetN(arg)!=1) {
    sprintf(errstr,"Expect scalar for %s",name);
    mexErrMsgTxt(errstr);
  }
  val=*mxGetPr(arg);
  if ((temp=floor(val))!=val) {
    sprintf(errstr,"Expect integer for %s",name);
    mexErrMsgTxt(errstr);
  }
  return (int) temp;
}

/*
 * This function also checks whether the element type is double. Use
 * 'getVecLenAnyType' to skip this check.
 */
int _getVecLen(const mxArray* arg,const char* name,char* errstr)
{
  int n;

  if (!mxIsDouble(arg)) {
    sprintf(errstr,"Expect real vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  if (mxGetM(arg)==0 || mxGetN(arg)==0)
    return 0;
  if ((n=mxGetM(arg))==1)
    n=mxGetN(arg);
  else if (mxGetN(arg)!=1) {
    sprintf(errstr,"Expect real vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  return n;
}

int _getVecLenAnyType(const mxArray* arg,const char* name,char* errstr)
{
  int n;

  if (mxGetM(arg)==0 || mxGetN(arg)==0)
    return 0;
  if ((n=mxGetM(arg))==1)
    n=mxGetN(arg);
  else if (mxGetN(arg)!=1) {
    sprintf(errstr,"Expect real vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  return n;
}

/*
 * NOTE: The string returned is allocated here using 'mxMalloc',
 * it has to be dealloc. by the user using 'mxFree'.
 */
const char* _getString(const mxArray* arg,const char* name,char* errstr)
{
  int len;
  char* buff;

  if (!mxIsChar(arg) || mxGetM(arg)!=1) {
    sprintf(errstr,"Expect char row vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  len=mxGetN(arg)+1;
  buff=(char*) mxMalloc(len*sizeof(char));
  mxGetString(arg,buff,len);

  return buff;
}

void _checkMatrix(const mxArray* arg,const char* name,int m,int n,char* errstr)
{
  if (!mxIsDouble(arg)) {
    sprintf(errstr,"Expect real matrix for %s",name);
    mexErrMsgTxt(errstr);
  }
  if (m!=-1 && mxGetM(arg)!=m) {
    sprintf(errstr,"Expect %d rows for %s",m,name);
    mexErrMsgTxt(errstr);
  }
  if (n!=-1 && mxGetN(arg)!=n) {
    sprintf(errstr,"Expect %d columns for %s",n,name);
    mexErrMsgTxt(errstr);
  }
}

void _checkVecIntNonneg(const mxArray* arg,const char* name,int n,char* errstr)
{
  int i,nn;
  double temp;
  const double* vP;

  if (!mxIsDouble(arg)) {
    sprintf(errstr,"Expect vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  if ((nn=mxGetM(arg))==1)
    nn=mxGetN(arg);
  else if (mxGetN(arg)!=1) {
    sprintf(errstr,"Expect vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  if (n!=-1 && nn!=n) {
    sprintf(errstr,"Expect length %d for %s",n,name);
    mexErrMsgTxt(errstr);
  }
  vP=mxGetPr(arg);
  for (i=0; i<nn; i++) {
    temp=*(vP++);
    if (floor(temp)!=temp || temp<0.0) {
      sprintf(errstr,"Expect nonnegative integer entries for %s",name);
      mexErrMsgTxt(errstr);
    }
  }
}

void _checkVecPosit(const mxArray* arg,const char* name,int n,char* errstr)
{
  int i,nn;
  double temp;
  const double* vP;

  if (!mxIsDouble(arg)) {
    sprintf(errstr,"Expect vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  if ((nn=mxGetM(arg))==1)
    nn=mxGetN(arg);
  else if (mxGetN(arg)!=1) {
    sprintf(errstr,"Expect vector for %s",name);
    mexErrMsgTxt(errstr);
  }
  if (n!=-1 && nn!=n) {
    sprintf(errstr,"Expect length %d for %s",n,name);
    mexErrMsgTxt(errstr);
  }
  vP=mxGetPr(arg);
  /*if (prefn!=-1 && prefn<nn) nn=prefn;*/
  for (i=0; i<nn; i++) {
    if (*(vP++)<=0.0) {
      sprintf(errstr,"Expect positive entries for %s",name);
      mexErrMsgTxt(errstr);
    }
  }
}

/*
 * Matrix argument can come with additional BLAS attributes. To pass these,
 * 'arg' must be a cell vector { BUFF, [YS XS M N], {SCODE} }. Here, BUFF
 * is a normal (buffer) matrix, YS, XS, M, N are integers, SCODE is a
 * string (optional). The matrix is BUFF(YS:(YS+M-1),XS:(XS+N-1)).
 *
 * Structure codes:
 * If SCODE is given, it must be a string of length 2. Pos.:
 * - 0: UPLO field, values 'U' (upper), 'L' (lower)
 * - 1: DIAG field, values 'U' (unit tri.), 'N' (non unit tri.)
 *      If UPLO spec., the def. value for DIAG is 'N'
 * A field value ' ' means: not specified.
 * These are passed to BLAS routines if required. The 'strcode' field
 * contains the codes separ. by 0, i.e. the C string for the codes
 * attached to each other.
 */
void _parseBLASMatrix(const mxArray* arg,const char* name,
		      fst_matrix* mat,int m,int n,char* errstr)
{
  int bm,bn,ys,xs,am,an,csz;
  const mxArray* bmat,*szvec,*scdvec;
  const double* iP;
  char sbuff[3];

  mat->strcode[0]=mat->strcode[2]=' ';
  mat->strcode[1]=mat->strcode[3]=0;
  if (!mxIsCell(arg)) {
    /* No cell array: Must be normal matrix */
    checkMatrix(arg,name,m,n);
    mat->buff=mxGetPr(arg);
    mat->stride=mat->m=mxGetM(arg); mat->n=mxGetN(arg);
  } else {
    if ((csz=mxGetM(arg)*mxGetN(arg))<2) {
      sprintf(errstr,"Array %s has wrong size",name);
      mexErrMsgTxt(errstr);
    }
    bmat=mxGetCell(arg,0);
    checkMatrix(bmat,name,-1,-1);
    bm=mxGetM(bmat); bn=mxGetN(bmat);
    if (getVecLen(mxGetCell(arg,1),name)!=4) {
      sprintf(errstr,"Index vector in %s has wrong size",name);
      mexErrMsgTxt(errstr);
    }
    iP=mxGetPr(mxGetCell(arg,1));
    ys=((int) iP[0])-1; xs=((int) iP[1])-1;
    am=(int) iP[2]; an=(int) iP[3];
    if (ys<0 || xs<0 || am<0 || an<0 || ys+am>bm || xs+an>bn) {
      sprintf(errstr,"Index vector in %s wrong",name);
      mexErrMsgTxt(errstr);
    }
    if ((m!=-1 && am!=m) || (n!=-1 && an!=n)) {
      sprintf(errstr,"Matrix %s has wrong size",name);
      mexErrMsgTxt(errstr);
    }
    mat->buff=mxGetPr(bmat)+(xs*bm+ys);
    mat->m=am; mat->n=an; mat->stride=bm;
    if (csz>2) {
      /* Structure codes */
      scdvec=mxGetCell(arg,2);
      if (!mxIsChar(scdvec) || mxGetM(scdvec)!=1 ||
	  mxGetN(scdvec)!=2) {
	sprintf(errstr,"Structure code string in %s wrong",name);
	mexErrMsgTxt(errstr);
      }
      mxGetString(scdvec,sbuff,3);
      if ((sbuff[0]!='U' && sbuff[0]!='L' && sbuff[0]!=' ') ||
	  (sbuff[1]!='U' && sbuff[1]!='N' && sbuff[1]!=' ')) {
	sprintf(errstr,"Structure code string in %s wrong",name);
	mexErrMsgTxt(errstr);
      }
      if (sbuff[0]!=' ' && sbuff[1]==' ')
	sbuff[1]='N'; /* def. value */
      if ((mat->m!=mat->n && sbuff[0]!=' ') ||
	  (sbuff[1]!=' ' && sbuff[0]==' ')) {
	sprintf(errstr,"Structure code string in %s inconsistent",name);
	mexErrMsgTxt(errstr);
      }
      mat->strcode[0]=sbuff[0]; mat->strcode[2]=sbuff[1];
    }
  }
}

/*
 * If 'n'!=-1, the length of the vector must be =='n'. If 'bn'!=-1, its
 * length must be <='bn'. In this case, if its size is <'bn', 'vec' represents
 * the initial part only.
 * NOTE: Striding value is set to 1.
 */
void _parseBLASVector(const mxArray* arg,const char* name,
		      fst_vector* vec,int n,int bn,char* errstr)
{
  int sz=getVecLen(arg,name);

  if (n!=-1 && sz!=n) {
    sprintf(errstr,"%s has wrong size",name);
    mexErrMsgTxt(errstr);
  }
  if (bn!=-1 && sz<bn) {
    sprintf(errstr,"%s is too short",name);
    mexErrMsgTxt(errstr);
    sz=(bn<n)?n:bn;
  }
  vec->n=sz; vec->buff=mxGetPr(arg);
  vec->stride=1;
}

/*
 * If 'n'!=-1, the length of the vector must be =='n'.
 * Length of vector is returned.
 */
int _parseDoubleVector(const mxArray* arg,const char* name,double** vec,int n,
		       char* errstr)
{
  int sz=getVecLen(arg,name);

  if (n!=-1 && sz!=n) {
    sprintf(errstr,"%s has wrong size",name);
    mexErrMsgTxt(errstr);
  }
  *vec=mxGetPr(arg);

  return sz;
}

/*
 * Vector must be of Matlab type int32, which corresponds to C++ type int.
 * If 'n'!=-1, the length of the vector must be =='n'.
 * Length of vector is returned.
 */
int _parseInt32Vector(const mxArray* arg,const char* name,int** vec,int n,
		      char* errstr)
{
  int sz=getVecLenAnyType(arg,name);

  if (!mxIsClass(arg,"int32")) {
    sprintf(errstr,"%s has wrong type (must be int32)",name);
    mexErrMsgTxt(errstr);
  }
  if (n!=-1 && sz!=n) {
    sprintf(errstr,"%s has wrong size",name);
    mexErrMsgTxt(errstr);
  }
  *vec=(int*) mxGetData(arg);

  return sz;
}

void** getZeroVoidArray(int n)
{
  int i;
  void** arr=(void**) mxMalloc(n*sizeof(void*));

  for (i=0; i<n; i++)
    arr[i]=(void*) 0;

  return arr;
}
