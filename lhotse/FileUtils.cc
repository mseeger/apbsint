//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Definition class FileUtils
 * ------------------------------------------------------------------- */

#include "lhotse/FileUtils.h"
#include "lhotse/NumberFormats.h"
#include <algorithm>

// Static members

const int FileUtils::boBigEndian;
const int FileUtils::boLittleEndian;
const int FileUtils::boOther;

bool FileUtils::tested(false);
int FileUtils::byteOrder(FileUtils::boOther);

void FileUtils::testFormats()
{
  if (!tested) {
    // Test assumptions
    if (sizeof(char)!=1)
      throw UnsupportedArchitectureException("Type char must have length 1 byte");
    int soi=sizeof(int);
    if (soi<4)
      throw UnsupportedArchitectureException("Type int must have length >=4 byte");
    cout << "FileUtils::testFormats:\n  Byte size int:     " << soi << endl
	 << "  Byte size double:  " << sizeof(double) << endl;

    // Test byte order on int
    int test=0x01020304; // might give compiler error on unsupp. architecture
    char* bt=(char*) &test;
    if (bt[soi-4+0]==1 && bt[soi-4+1]==2 && bt[soi-4+2]==3 && bt[soi-4+3]==4) {
      byteOrder=boBigEndian;
      cout << "  System bute order: big endian." << endl;
    } else if (bt[0]==4 && bt[1]==3 && bt[2]==2 && bt[3]==1) {
      byteOrder=boLittleEndian;
      cout << "  System byte order: little endian." << endl;
    } else {
      byteOrder=boOther;
      cout << "  WARNING: System byte order neither big nor little endian!"
	   << endl;
    }
    tested=true;
  }
}

void FileUtils::saveSeq(ofstream& os,const char* data,int size,int n,int step,
			const int* index,int fsize)
{
  int i,j;
  ArrayHandle<char> buff;

  if (size<1) throw InvalidParameterException(EXCEPT_MSG("size"));
  if (fsize==-1) fsize=size;
  else if (fsize<1) throw InvalidParameterException(EXCEPT_MSG("Fsize"));
  testFormats();
  bool revOrder=(size>1 && byteOrder==boLittleEndian);
  if (n>0) {
    if (index==0) {
      if (step==0) throw InvalidParameterException(EXCEPT_MSG("step"));
      if (fsize==size && !revOrder) {
	if (step==1)
	  os.write(data,n*size);
	else {
	  for (i=0; i<n; i++,data+=size*step)
	    os.write(data,size);
	}
      } else if (fsize>=size) {
	buff.changeRep(fsize);
	std::fill(buff.p(),buff.p()+fsize,0);
	for (i=0; i<n; i++,data+=size*step) {
	  if (revOrder)
	    reverseBO(buff.p()+(fsize-size),data,size);
	  else
	    memmove(buff.p()+(fsize-size),data,size);
	  os.write(buff.p(),fsize);
	}
      } else {
	buff.changeRep(size);
	for (i=0; i<n; i++,data+=size*step) {
	  if (revOrder) reverseBO(buff.p(),data,size);
	  else memmove(buff.p(),data,size);
	  for (j=0; j<size-fsize; j++)
	    if (buff[j]!=0)
	      throw FileFormatException(EXCEPT_MSG("Byte size insufficient"));
	  os.write(buff.p()+(size-fsize),fsize);
	}
      }
    } else {
      if (fsize==size && !revOrder) {
	for (i=0; i<n; i++)
	  os.write(data+((*(index++))*size),size);
      } else if (fsize>=size) {
	buff.changeRep(fsize);
	std::fill(buff.p(),buff.p()+fsize,0);
	for (i=0; i<n; i++) {
	  if (revOrder)
	    reverseBO(buff.p()+(fsize-size),data+((*(index++))*size),size);
	  else
	    memmove(buff.p()+(fsize-size),data+((*(index++))*size),size);
	  os.write(buff.p(),fsize);
	}
      } else {
	buff.changeRep(size);
	for (i=0; i<n; i++) {
	  if (revOrder) reverseBO(buff.p(),data+((*(index++))*size),size);
	  else memmove(buff.p(),data+((*(index++))*size),size);
	  for (j=0; j<size-fsize; j++)
	    if (buff[j]!=0)
	      throw FileFormatException(EXCEPT_MSG("Byte size insufficient"));
	  os.write(buff.p()+(size-fsize),fsize);
	}
      }
    }
  }
}

void FileUtils::loadSeq(ifstream& is,char* data,int size,int n,int step,
			const int* index,int fsize)
{
  int i,j;
  ArrayHandle<char> buff;

  if (size<1) throw InvalidParameterException(EXCEPT_MSG("size"));
  if (fsize==-1) fsize=size;
  else if (fsize<1)
    throw InvalidParameterException(EXCEPT_MSG("fsize"));
  testFormats();
  bool revOrder=(size>1 && byteOrder==boLittleEndian);
  if (n>0) {
    if (index==0) {
      if (step==0)
	throw InvalidParameterException(EXCEPT_MSG("step"));
      if (size==fsize && !revOrder) {
	if (step==1)
	  is.read(data,n*size);
	else {
	  for (i=0; i<n; i++,data+=step*size)
	    is.read(data,size);
	}
      } else if (fsize<=size) {
	buff.changeRep(size);
	std::fill(buff.p(),buff.p()+size,0);
	for (i=0; i<n; i++,data+=size*step) {
	  // File is always big endian, so leading bytes 0
	  is.read(buff.p()+(size-fsize),fsize);
	  if (revOrder) reverseBO(data,buff.p(),size);
	  else memmove(data,buff.p(),size);
	}
      } else {
	buff.changeRep(fsize);
	for (i=0; i<n; i++,data+=size*step) {
	  // File is always big endian. Skip msb's, check whether 0
	  is.read(buff.p(),fsize);
	  for (j=0; j<fsize-size; j++)
	    if (buff[j]!=0)
	      throw FileFormatException(EXCEPT_MSG("Byte size insufficient"));
	  if (revOrder) reverseBO(data,buff.p()+(fsize-size),size);
	  else memmove(data,buff.p()+(fsize-size),size);
	}
      }
    } else {
      if (size==fsize && !revOrder) {
	for (i=0; i<n; i++)
	  is.read(data+((*(index++))*size),size);
      } else if (fsize<=size) {
	ArrayHandle<char> buff(size);
	std::fill(buff.p(),buff.p()+size,0);
	for (i=0; i<n; i++) {
	  // File is always big endian, so leading bytes 0
	  is.read(buff.p()+(size-fsize),fsize);
	  if (revOrder)
	    reverseBO(data+((*(index++))*size),buff.p(),size);
	  else
	    memmove(data+((*(index++))*size),buff.p(),size);
	}
      } else {
	buff.changeRep(fsize);
	for (i=0; i<n; i++,data+=size*step) {
	  // File is always big endian. Skip msb's. Check whether 0
	  is.read(buff.p(),fsize);
	  for (j=0; j<fsize-size; j++)
	    if (buff[j]!=0)
	      throw FileFormatException(EXCEPT_MSG("Byte size insufficient"));
	  if (revOrder)
	    reverseBO(data+((*(index++))*size),buff.p()+(fsize-size),size);
	  else
	    memmove(data+((*(index++))*size),buff.p()+(fsize-size),size);
	}
      }
    }
  }
}

void FileUtils::saveBoolCompact(ofstream& os,const bool* data,int size)
{
  int i,mod=0;
  uchar byte=0;

  if (size<0) throw InvalidParameterException(EXCEPT_MSG(""));
  for (i=0; i<size; i++) {
    byte<<=1;
    if (data[i]) byte&=1;
    if (++mod==8) {
      os.write((char*) &byte,1);
      byte=0; mod=0;
    }
  }
  if (mod!=0) os.write((char*) &byte,1);
}

void FileUtils::loadBoolCompact(ifstream& is,bool* data,int size)
{
  int i,limit,j;
  uchar byte,mask;

  if (size<0) throw InvalidParameterException(EXCEPT_MSG(""));
  limit=size/8;
  for (i=j=0; i<limit; i++) {
    is.read((char*) &byte,1);
    for (mask=0x80; mask!=0; mask>>=1)
      data[j++]=((byte&mask)!=0);
  }
  if (j<size) {
    is.read((char*) &byte,1);
    for (mask=0x80>>(8-size+j); j<size; j++,mask>>=1)
      data[j]=((byte&mask)!=0);
  }
}

int FileUtils::loadHeader(ifstream& is,const string& tag,bool addAdd,
			  bool noVer)
{
  int len=tag.length()+(addAdd?1:0);
  char buff[len+1];

  is.read(buff,len); buff[len]=0;
  string tempTag(addAdd?"@":""); tempTag+=tag;
  if (strcmp(buff,tempTag.c_str())!=0) {
    string msg="Unknown tag. Expected: '";
    msg+=tempTag; msg+="'";
    throw FileFormatException(msg.c_str());
  }
  int ffVer;
  if (!noVer)
    NumberFormats<int>::load(is,&ffVer,1,1,0,4);
  else
    noVer=0;

  return ffVer;
}

int FileUtils::loadHeaderFlex(ifstream& is,const string& tag,bool* oldFormat,
			      bool noVer)
{
  int len=tag.length(),actlen;
  char buff[len];
  char* act;

  is.read(buff,1);
  if (buff[0]=='@') {
    actlen=len;
    act=buff;
  } else {
    actlen=len-1;
    act=buff+1;
  }
  is.read(act,actlen); buff[len]=0;
  if (strcmp(buff,tag.c_str())!=0) {
    string msg="Unknown tag. Expected: '";
    msg+=tag+"' or '@";
    msg+=tag+"'";
    throw FileFormatException(msg.c_str());
  }
  int ffVer;
  if (!noVer)
    NumberFormats<int>::load(is,&ffVer,1,1,0,4);
  else
    noVer=0;

  return ffVer;
}

int FileUtils::loadHeaderMulti(ifstream& is,const ArrayHandle<string>& tagList,
			       int& resInd,bool noVer)
{
  int i,numOK=tagList.size(),pos=0;
  ArrayHandle<int> okPos(numOK);
  char act;

  if (numOK==0) throw InvalidParameterException("tagList");
  for (i=0; i<numOK; i++) okPos[i]=i;
  while (numOK>1) {
    is.read(&act,1);
    for (i=0; i<numOK; )
      if (act!=tagList[okPos[i]][pos]) {
	okPos[i]=okPos[numOK-1];
	numOK--;
      } else i++;
    pos++;
  }
  if (numOK==0) throw FileFormatException("Unknown file tag");
  int cand=okPos[0];
  int len=tagList[cand].length();
  if (len>pos) {
    ArrayHandle<char> buff(len-pos+1);
    is.read(buff.p(),len-pos); buff[len-pos]=0;
    if (strcmp(buff.p(),tagList[cand].substr(pos).c_str())!=0)
      throw FileFormatException("Unknown file tag");
  }
  int ffVer;
  if (!noVer)
    NumberFormats<int>::load(is,&ffVer,1,1,0,4);
  else
    noVer=0;

  resInd=cand;
  return ffVer;
}

void FileUtils::saveHeader(ofstream& os,const string& tag,int ffVer,
			   bool addAdd)
{
  if (addAdd)
    os.write("@",1);
  os.write(tag.c_str(),tag.length());
  NumberFormats<int>::save(os,&ffVer,1,1,0,4);
}
