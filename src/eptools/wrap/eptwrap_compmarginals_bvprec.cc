/* -------------------------------------------------------------------
 * EPTWRAP_COMPMARGINALS_BVPREC
 *
 * Computes marginals of a and c parameters for models with bivariate
 * precision potentials.
 * Message parameters are RP_A, RP_C, marginals written to MARGA,
 * MARGC. The first size(RP_A) entries of RP_TAUIND is the index
 * [k(j)], mapping RP_{A|C} to MARG{A|C} entries.
 * NOTE: RP_TAUIND can have further entries, which are not used here.
 *
 * Input:
 * - RP_TAUIND Index k(j) [int32 array]
 * - RP_A:     Message parameters [double array]
 * - RP_C:     "
 * - MARGA:    Marginal parameters written here [double array]
 * - MARGC:    "
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include "src/main.h"
#include "src/eptools/wrap/eptools_helper.h"
#include "src/eptools/wrap/eptwrap_compmarginals_bvprec.h"
#include <algorithm>

void eptwrap_compmarginals_bvprec(int ain,int aout,W_IARRAY(rp_tauind),
				  W_DARRAY(rp_a),W_DARRAY(rp_c),
				  W_DARRAY(marga),W_DARRAY(margc),W_ERRORARGS)
{
  try {
    /* Read arguments */
    if (ain!=5)
      W_RETERROR(2,"Need 5 input arguments");
    if (aout!=0)
      W_RETERROR(2,"No return arguments");
    W_CHKSIZE(rp_c,nrp_a,"RP_C");
    if (nrp_tauind < nrp_a)
      W_RETERROR(1,"RP_TAUIND shorter than RP_A");
    W_CHKSIZE(margc,nmarga,"MARGC");
    std::fill(marga,marga+nmarga,0.0);
    std::fill(margc,margc+nmarga,0.0);
    int i,j;
    for (i=0; i<nrp_a; i++) {
      marga[j=rp_tauind[i]] += rp_a[i];
      margc[j] += rp_c[i];
    }
    W_RETOK;
  } catch (StandardException ex) {
    W_RETERROR_ARGS(1,"Caught LHOTSE exception: %s",ex.msg());
  } catch (...) {
    W_RETERROR(1,"Caught unspecified exception");
  }
}
