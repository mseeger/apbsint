/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class SpecfunServices
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_SPECFUNSERVICES_H
#define EPTOOLS_SPECFUNSERVICES_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/default.h"

//BEGINNS(eptools)
  /**
   * Collects static methods for computing certain special functions.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
#ifndef HAVE_LIBGSL
#include "src/eptools/potentials/SpecfunServices_basic.h"
#else
#include "src/eptools/potentials/SpecfunServices_workaround.h"
#endif
//ENDNS

#endif
