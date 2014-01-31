/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Definition of class EPPotentialNamedFactory
 * ------------------------------------------------------------------- */

#include "src/eptools/potentials/EPPotentialNamedFactory.h"

//BEGINNS(eptools)
  MAP_TYPE(string,int) EPPotentialNamedFactory::potNames;
  MAP_TYPE(int,string) EPPotentialNamedFactory::potIDs;

  void EPPotentialNamedFactory::setup()
  {
    if (potNames.empty()) {
      potNames["Gaussian"]     = potGaussian;
      potIDs[potGaussian]      = "Gaussian";
      potNames["Laplace"]      = potLaplace;
      potIDs[potLaplace]       = "Laplace";
      potNames["Probit"]       = potProbit;
      potIDs[potProbit]        = "Probit";
      potNames["Heaviside"]    = potHeaviside;
      potIDs[potHeaviside]     = "Heaviside";
      potNames["Exponential"]  = potExponential;
      potIDs[potExponential]   = "Exponential";
      potNames["QuantRegress"] = potQuantRegress;
      potIDs[potQuantRegress]  = "QuantRegress";
      potNames["GaussMixture"] = potGaussMixture;
      potIDs[potGaussMixture]  = "GaussMixture";
      potNames["SpikeSlab"]    = potSpikeSlab;
      potIDs[potSpikeSlab]     = "SpikeSlab";
    }
  }
//ENDNS
