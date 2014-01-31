/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotentialNamedFactory
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTENTIALNAMEDFACTORY_H
#define EPTOOLS_EPPOTENTIALNAMEDFACTORY_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPPotentialFactory.h"

//BEGINNS(eptools)
  /**
   * Extends 'EPPotentialFactory' by association
   *
   *   Name (string) -> ID (nonneg. int)
   *
   * External interfaces refer to EP potentials by Name, but
   * internal representations typically translate this to ID.
   * This functionality is not in 'EPPotentialFactory', which should
   * be as lean as possible (is created on the fly for Matlab MEX
   * calls).
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotentialNamedFactory : public EPPotentialFactory
  {
  protected:
    // Additional members

    static MAP_TYPE(string,int) potNames;
    static MAP_TYPE(int,string) potIDs;

  public:
    // Public static methods

    /**
     * @param name Name string
     * @return     ID, or -1 if not in 'potNames'
     */
    static int getID4Name(const string& name) {
      setup();
      MAP_CONSTITER(string,int) it=potNames.find(name);
      if (it!=potNames.end())
	return it->second;
      else
	return -1;
    }

    /**
     * Throws 'OutOfRangeException' if 'pid' is not a valid ID.
     *
     * @param pid ID number
     * @return    Potential name
     */
    static const string& getName4ID(int pid) {
      setup();
      if (!isValidID(pid)) throw OutOfRangeException(EXCEPT_MSG(""));

      return potIDs[pid];
    }

    static EPScalarPotential* create(int pid,const double* pv) {
      return EPPotentialFactory::create(pid,pv);
    }

    static EPScalarPotential* create(const string& name,const double* pv) {
      setup();
      int pid=getID4Name(name);
      if (pid==-1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      return EPPotentialFactory::create(pid,pv);
    }

    static EPScalarPotential* createDefault(int pid,const double* pv=0) {
      return EPPotentialFactory::createDefault(pid,pv);
    }

    static EPScalarPotential* createDefault(const string& name,
					    const double* pv=0) {
      setup();
      int pid=getID4Name(name);
      if (pid==-1)
	throw InvalidParameterException(EXCEPT_MSG(""));
      return EPPotentialFactory::createDefault(pid,pv);
    }

  protected:
    /**
     * If 'potNames' is empty, it is setup. Otherwise, do nothing. This
     * method is called by all others.
     */
    static void setup();
  };
//ENDNS

#endif
