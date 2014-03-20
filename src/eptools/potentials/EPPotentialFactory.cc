/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotentialFactory
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/EPPotentialFactory.h"
#include "src/eptools/potentials/EPPotGaussian.h"
#include "src/eptools/potentials/EPPotProbit.h"
#include "src/eptools/potentials/EPPotQuantileRegress.h"
#include "src/eptools/potentials/EPPotLaplace.h"
#include "src/eptools/potentials/EPPotGaussMixture.h"
#include "src/eptools/potentials/EPPotSpikeSlab.h"

//BEGINNS(eptools)
  const int EPPotentialFactory::potGaussian;
  const int EPPotentialFactory::potLaplace;
  const int EPPotentialFactory::potProbit;
  const int EPPotentialFactory::potHeaviside;
  const int EPPotentialFactory::potExponential;
  const int EPPotentialFactory::potQuantRegress;
  const int EPPotentialFactory::potGaussMixture;
  const int EPPotentialFactory::potSpikeSlab;
  const int EPPotentialFactory::potLast;
#ifdef HAVE_WORKAROUND
#include "src/eptools/potentials/EPPotentialFactory_workaround.cc"
#endif

  /*
   * We call the static method 'getArgumentGroup_static'.
  */
  int EPPotentialFactory::getArgumentGroup(int pid)
  {
    int ret=0;

    if (!isValidID(pid))
      throw InvalidParameterException(EXCEPT_MSG(""));
    switch (pid) {
    case potGaussian:
      ret = EPPotGaussian::getArgumentGroup_static();
      break;
    case potLaplace:
      ret = EPPotLaplace::getArgumentGroup_static();
      break;
    case potProbit:
    case potHeaviside:
      ret = EPPotProbit::getArgumentGroup_static();
      break;
    case potExponential:
      throw NotImplemException(EXCEPT_MSG(""));
      break;
    case potQuantRegress:
      ret = EPPotQuantileRegress::getArgumentGroup_static();
      break;
    case potGaussMixture:
      ret = EPPotGaussMixture::getArgumentGroup_static();
      break;
    case potSpikeSlab:
      ret = EPPotSpikeSlab::getArgumentGroup_static();
      break;
#ifdef HAVE_WORKAROUND
    default:
      ret = getArgumentGroup_workaround(pid);
#endif
    }

    return ret;
  }

  EPScalarPotential* EPPotentialFactory::create(int pid,const double* pv,
						void* annot)
  {
    int i;
    EPScalarPotential* rpot=0;
    EPPotGaussMixture* potGM=0;

    if (!isValidID(pid) || pv==0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    switch (pid) {
    case potGaussian:
      rpot = new EPPotGaussian(pv[0],pv[1]);
      break;
    case potLaplace:
      rpot = new EPPotLaplace(pv[0],pv[1]);
      break;
    case potProbit:
      rpot = new EPPotProbit(pv[0],pv[1],false);
      break;
    case potHeaviside:
      rpot = new EPPotProbit(pv[0],pv[1],true);
      break;
    case potExponential:
      throw NotImplemException(EXCEPT_MSG(""));
      break;
    case potQuantRegress:
      rpot = new EPPotQuantileRegress(pv[0],pv[1],pv[2]);
      break;
    case potGaussMixture:
      i=(int) ceil(pv[0]);
      rpot = potGM = new EPPotGaussMixture(i);
      potGM->setCVals(pv+1);
      potGM->setVariances(pv+i);
      break;
    case potSpikeSlab:
      rpot = new EPPotSpikeSlab(pv[0],pv[1]);
      break;
#ifdef HAVE_WORKAROUND
    default:
      rpot = create_workaround(pid,pv,annot);
#endif
    }

    return rpot;
  }

  EPScalarPotential* EPPotentialFactory::createDefault(int pid,const double* pv,
						       void* annot)
  {
    int i;
    EPScalarPotential* rpot=0;

    if (!isValidID(pid) || pv==0)
      throw InvalidParameterException(EXCEPT_MSG(""));
    switch (pid) {
    case potGaussian:
      rpot = new EPPotGaussian();
      break;
    case potLaplace:
      rpot = new EPPotLaplace();
      break;
    case potProbit:
      rpot = new EPPotProbit(false);
      break;
    case potHeaviside:
      rpot = new EPPotProbit(true);
      break;
    case potExponential:
      throw NotImplemException(EXCEPT_MSG(""));
      break;
    case potQuantRegress:
      rpot = new EPPotQuantileRegress();
      break;
    case potGaussMixture:
      if (pv==0)
	throw InvalidParameterException("Need construction parameters");
      i=(int) ceil(pv[0]);
      rpot = new EPPotGaussMixture(i);
      break;
    case potSpikeSlab:
      rpot = new EPPotSpikeSlab();
      break;
#ifdef HAVE_WORKAROUND
    default:
      rpot = createDefault_workaround(pid,pv,annot);
#endif
    }

    return rpot;
  }
//ENDNS
