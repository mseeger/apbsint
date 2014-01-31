//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Header class FileUtils
 * ------------------------------------------------------------------- */

/*
 * NOTE: Not a completely portable solution! Files may have to be converted
 * explicitly for architectures not supported here.
 * ==> TODO: Use serialization library to do this!
 */

#ifndef FILEUTILS_H
#define FILEUTILS_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "lhotse/global.h" // global header

template<class T> class NumberFormats;

/**
 * Provides static methods for basic file access.
 * The characteristics of the basic number formats of the architecture are
 * tested automatically (char,int,double). Methods are provided to store
 * objects in normal or reversed byte order.
 * <p>
 * ATTENTION: These methods are only tested on SOME architectures. When
 * porting LHOTSE to a yet unsupported architecture, it is required to
 * test the number formats there, and to check whether the methods here
 * are appropriate or have to be modified.
 * <p>
 * There is a source architecture (file written) and a target architecture
 * (file read). Assumptions:
 * 1. char has size of 1 byte
 * 2. int has at least 4 bytes
 * 3. Source and target architecture have the same byte sizes for all
 *    types
 *    ==> Major source of non-portability!
 * 4. Each architecture might have a certain byte order, i.e. the ordering
 *    in which multi-byte types are stored in memory. This can be:
 *    - big endian: most signif. byte at lowest address
 *    - little endian: least signif. byte at lowest address
 *    - another (not really supported)
 *    We assume that the SAME byte order is applied to ANY multi-byte type
 *    (here: int,float,double)
 * 5. For a byte order other than big or little endian, we just use the
 *    system byte order
 *    ==> Source of non-portability
 * 6. The format of the multi-byte types vary across supported architectures
 *    ONLY BY THEIR BYTE ORDER
 * The assumptions are checked in 'testFormats'. If they are not
 * met, an exception is thrown.
 * NOTE: It is STRONGLY discouraged to use LHOTSE in the present form on
 * an architecture with a byte order different from big/little endian.
 * The file formats typically do NOT indicate whether they have been written
 * in a machine-dep. way because of this, and attempts to load them on other
 * architectures result in undetected errors.
 * Instead, this class should be ported to the architecture by introducing a
 * new byte order code and modifying ALL methods.
 * <p>
 * To achieve portability amongst the supported architectures, we generally
 * use the big endian byte order in files. This is the favourite byte order
 * for UNIX systems. In the moment, an architecture is called "supported"
 * w.r.t. this class if it satisfies the assumptions above and has big or
 * little endian byte order, the same for all multi-byte types.
 * The byte order is tested on int only.
 * <p>
 * Different byte sizes:
 * Types like int may have a different byte size on different systems. For
 * portability, the byte size of used types must be fixed for a file format,
 * or it must be stored in the file header.
 * The byte size to be used in files can be passed to 'loadSeq', 'saveSeq'
 * in the 'fsize' arg. This is relevant for the corr. 'NumberFormats'
 * wrappers, where the system byte size is passed in 'size' (det. by
 * sizeof(T)), which can be different from 'fsize'.
 * NOTE: If 'fsize' < 'size' in 'saveSeq' or 'size' < 'fsize' in 'loadSeq',
 * the msb's are chopped. The methods check whether the chopped msb's are
 * 0, and throw an exception if not.
 * <p>
 * NOTE: The methods for storing/loading arrays def. here should not be used
 * directly, use the 'NumberFormats' methods instead.
 *
 * @author  Matthias Seeger
 * @version %I% %G%
 */
class FileUtils
{
  template<class T> friend class NumberFormats;

protected:
  // Static members, constants

  static bool tested;   // has architecture already been tested?
  static int byteOrder; // byte order

  static const int boBigEndian   =0; // Byte orders
  static const int boLittleEndian=1;
  static const int boOther       =2;

public:
  // Public static methods

  /**
   * Tests number formats on the current architecture and sets 'byteOrder'.
   * This is done only if 'tested' is still false. The method sets 'tested'
   * to true, therefore the tests are never done twice in a program run.
   * If any assumptions are violated, we throw a
   * 'UnsupportedArchitectureException' exception.
   */
  static void testFormats();

  /**
   * Returns byte order
   *
   * @return Byte order (use constants 'boXXX' def. above)
   */
  static int getByteOrder() {
    testFormats();
    return byteOrder;
  }

  /**
   * Stores sequence of bool in a machine-independent way, compact way.
   * If b_1,...,b_n is the sequence, we store bytes
   *   (b_1 ... b_8), ..., (0 ... 0 b_{8k+1} .. b_n),
   * where k = floor((n-1)/8).
   *
   * @param os   Output file stream
   * @param data Pointer to bool sequence
   * @param size Length of sequence
   */
  static void saveBoolCompact(ofstream& os,const bool* data,int size);

  /**
   * Loads sequence of bool, as written by 'saveBoolCompact'.
   *
   * @param is   Input file stream
   * @param data Pointer to buffer for bool sequence
   * @param size Length of sequence
   */
  static void loadBoolCompact(ifstream& is,bool* data,int size);

  /**
  * Opens file with name 'fname' for reading, using input filestream variable
  * 'is'. In case of an error, a 'FileUtilsException' is thrown containing the
  * name of the file.
  *
  * @param fname Filename
  * @param is    Input filestream
  */
  static void openFileRead(const char* fname,ifstream& is) {
    is.open(fname);
    if (!is) {
      string msg="Cannot open file '";
      msg+=fname; msg+="' for reading!";
      throw FileUtilsException(msg.c_str());
    }
  }

  static void openFileRead(const string& fname,ifstream& is) {
    openFileRead(fname.c_str(),is);
  }

  /**
   * Creates file with name 'fname' for writing, using output filestream
   * variable 'os'. In case of an error, a 'FileUtilsException' is thrown
   * containing the name of the file.
   *
   * @param fname Filename
   * @param os    Output filestream
   */
  static void openFileWrite(const char* fname,ofstream& os) {
    os.open(fname);
    if (!os) {
      string msg="Cannot create file '";
      msg+=fname; msg+="' for writing!";
      throw FileUtilsException(msg.c_str());
    }
  }

  static void openFileWrite(const string& fname,ofstream& os) {
    openFileWrite(fname.c_str(),os);
  }

  /**
   * Opens file with name 'fname' for appending, using stream 'os'.
   * In case of an error, a 'FileUtilsException' is thrown
   * containing the name of the file.
   *
   * @param fname Filename
   * @param os    Output filestream
   */
  static void openFileAppend(const char* fname,ofstream& os) {
    os.open(fname,ios::app);
    if (!os) {
      string msg="Cannot open file '";
      msg+=fname; msg+="' for appending!";
      throw FileUtilsException(msg.c_str());
    }
  }

  static void openFileAppend(const string& fname,ofstream& os) {
    openFileAppend(fname.c_str(),os);
  }

  /**
   * Loads header of LHOTSE object file (stream 'is'), cons. of tag and FF
   * version number. This method applies to binary object files using the
   * current LHOTSE convention:
   * - tag [sequence of char, no 0-termination]
   * - FF version number [int]
   * Usually, tags have the form '@<classname>'. The tag must be passed in
   * 'tag'. The methods first reads a char seq. of the length of 'tag' and
   * compares it against 'tag'. If they are not equal, a 'FileFormatException'
   * is thrown. Otherwise, the method reads one int (the FF version number)
   * and returns it.
   * If 'addAdd' is true (def.: false), a leading '@' is added: we read
   * a char seq. of length('tag')+1, compare the first char. against '@' and
   * the rest against 'tag'.
   * <p>
   * NOTE: Some older LHOTSE file formats have different header formats,
   * this method cannot be used on them! See also 'loadHeaderFlex'.
   * <p>
   * If 'noVer' is true, the FF version number is not read, and 0 is returned.
   *
   * @param is     Input file stream
   * @param tag    File tag (to be compared with true file tag)
   * @param addAdd Optional. Def.: false. See above
   * @param noVer  Optional. Def.: false. See above
   * @return       FF version number
   */
  static int loadHeader(ifstream& is,const string& tag,bool addAdd=false,
			bool noVer=false);

  /**
   * Loads header of LHOTSE object file (stream 'is'), cons. of tag and FF
   * version number. This method applies to binary object files using the
   * current LHOTSE convention:
   * - tag [sequence of char, no 0-termination]
   * - FF version number [int]
   * By LHOTSE convention, tags have the form '@<classname>', however some
   * older formats have '<classname>' instead. This method should be used on
   * classes for which such older formats coexist with newer formats. We
   * first read the tag and compare it to 'tag' and to the concat. of '@' with
   * 'tag'. If none of these match, a 'FileFormatException' is thrown.
   * Otherwise, the method reads one int (the FF version number) and returns
   * it. Optionally, in 'oldFormat', we return true iff the tag is equal to
   * 'tag' without the '@' prefix.
   * <p>
   * If 'noVer' is true, the FF version number is not read, and 0 is returned.
   *
   * @param is        Input file stream
   * @param tag       File tag. See above
   * @param oldFormat Optional. See above
   * @param noVer     See above. Def.: false
   * @return          FF version number
   */
  static int loadHeaderFlex(ifstream& is,const string& tag,bool* oldFormat=0,
			    bool noVer=false);

  /**
   * Loads header of LHOTSE object file (stream 'is'), cons. of tag and FF
   * version number. This method applies to binary object files using the
   * current LHOTSE convention:
   * - tag [sequence of char, no 0-termination]
   * - FF version number [int]
   * A list of potential tags is passed in 'tagList'. None of the entries must
   * be the prefix of another one. The method reads characters until only one
   * of the entries matches. If none of the entries match, a
   * 'FileFormatException' is thrown.
   * Otherwise, the method reads one int (the FF version number) and returns
   * it. The index of the file tag within 'tagList' is ret. via 'resInd'.
   * <p>
   * If 'noVer' is true, the FF version number is not read, and 0 is returned.
   *
   * @param is      Input file stream
   * @param tagList List of pot. file tag (to be compared with true file tag)
   * @param resInd  See above
   * @param noVer   See above. Def.: false
   * @return        FF version number
   */
  static int loadHeaderMulti(ifstream& is,const ArrayHandle<string>& tagList,
			     int& resInd,bool noVer=false);	

  /**
   * Stores header for a LHOTSE object file (into stream 'os'). See
   * 'loadHeader' for the format. The tag must be passed in 'tag', the FF
   * version number in 'ffVer'.
   * If 'addAdd' is true (def.: false), a leading '@' is added to 'tag' before
   * writing.
   *
   * @param os     Output file stream
   * @param tag    File tag (typically @<classname>)
   * @param ffVer  FF version number
   * @param addAdd Optional. Def.: false. See above
   */
  static void saveHeader(ofstream& os,const string& tag,int ffVer,
			 bool addAdd=false);

protected:
  // Internal static methods (inline)

  /**
   * Stores sequence of length 'n', starting from 'data'. The sequence consists
   * of elements of 'size' bytes each. For byte order little endian, the order
   * of bytes is reversed for each element before the byte stream is written.
   * If 'index' is given, then 'n' is the length of the index, and 'data' is
   * indexed in the sense that element i is found at
   *   data[index[i]*size],...,data[index[i]*(size+1)-1].
   * If 'index' is not given, a step size can be passed via 'step' (!=0),
   * coding the index i*step. 'step' can be negative.
   * <p>
   * Different byte size in file stream:
   * If 'fsize' is given, the entries in the file stream 'is' have byte
   * size 'fsize', which can be different from 'size'. See header comment.
   *
   * @param os    Output file stream
   * @param data  Pointer to sequence
   * @param size  Size of elements in byte
   * @param n     Length of sequence. Def.: 1
   * @param step  Step size. Ignored if 'index' given. Def.: 1
   * @param index See above. Def.: 0
   * @param fsize S.a. Def.: -1 (-> 'size')
   */
  static void saveSeq(ofstream& os,const char* data,int size,int n=1,
		      int step=1,const int* index=0,int fsize=-1);

  /**
   * Loads sequence of length 'n', writing into 'data'. The sequence consists
   * of elements of 'size' bytes each. For byte order little endian, the order
   * of bytes is reversed for each element before 'data' is written.
   * If 'index' is given, then 'n' is the length of the index, and 'data' is
   * indexed in the sense that element i is written to
   *   data[index[i]*size],...,data[index[i]*(size+1)-1]
   * If 'index' is not given, a step size can be passed via 'step' (!=0),
   * coding the index i*step. 'step' can be negative.
   * <p>
   * Different byte size in file stream:
   * If 'fsize' is given, the entries in the file stream 'is' have byte
   * size 'fsize', which can be different from 'size'. See header comment.
   *
   * @param is    Input file stream
   * @param data  Pointer to target buffer
   * @param size  Size of elements (in byte) in 'data'
   * @param n     Length of sequence. Def.: 1
   * @param step  Step size. Ignored if 'index' given. Def.: 1
   * @param index See above. Def.: 0
   * @param fsize S.a. Def.: -1 (-> 'size')
   */
  static void loadSeq(ifstream& is,char* data,int size,int n=1,int step=1,
		      const int* index=0,int fsize=-1);

  /**
   * Reverses order of byte sequence 'src', writes into 'trg'. The buffers
   * must be different!
   *
   * @param trg Target buffer
   * @param src Source buffer
   */
  static void reverseBO(char* trg,const char* src,int size) {
    switch(size) {
    case 8:
      trg[0]=src[7]; trg[1]=src[6]; trg[2]=src[5]; trg[3]=src[4];
      trg[4]=src[3]; trg[5]=src[2]; trg[6]=src[1]; trg[7]=src[0];
      break;
    case 4:
      trg[0]=src[3]; trg[1]=src[2]; trg[2]=src[1]; trg[3]=src[0];
      break;
    case 2:
      trg[0]=src[1]; trg[1]=src[0];
      break;
    default:
      register int i=0;
      const char* src2=src+(size-1);
      for (i=0; i<size; i++) *(trg++)=*(src2--);
      break;
    }
  }
};

#endif
