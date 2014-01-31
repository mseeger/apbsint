//NEWCODE
/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Library source file
 * Module: GLOBAL
 * Desc.:  Definition class Interval
 * ------------------------------------------------------------------- */

#include "lhotse/global.h"
#include "lhotse/Interval.h"

template<> bool DefIVal<char>::isInit(false);
template<> bool DefIVal<uchar>::isInit(false);
template<> bool DefIVal<int>::isInit(false);
template<> bool DefIVal<uint>::isInit(false);
template<> bool DefIVal<long>::isInit(false);
template<> bool DefIVal<unsigned long>::isInit(false);
template<> bool DefIVal<float>::isInit(false);
template<> bool DefIVal<double>::isInit(false);

template<> char DefIVal<char>::zeroVal(0);
template<> uchar DefIVal<uchar>::zeroVal(0);
template<> int DefIVal<int>::zeroVal(0);
template<> uint DefIVal<uint>::zeroVal(0);
template<> long DefIVal<long>::zeroVal(0);
template<> unsigned long DefIVal<unsigned long>::zeroVal(0);
template<> float DefIVal<float>::zeroVal(0.0);
template<> double DefIVal<double>::zeroVal(0.0);
