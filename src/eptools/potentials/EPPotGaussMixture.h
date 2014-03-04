/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotGaussMixture
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTGAUSSMIXTURE_H
#define EPTOOLS_EPPOTGAUSSMIXTURE_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"
#include <algorithm>

//BEGINNS(eptools)
  /**
   * Mixture of Gaussians potential:
   *   t(s) = sum_{l=0}^{L-1} p_l N(s|0,v_l),   L>=2, v_l>0.
   * Here:
   *   p_l = exp(c_l) / [sum_k exp(c_k)]   and c_{L-1} = 0
   * Parameters: L, c_0, ..., c_{L-2}, v_0, ..., v_{L-1} [2*L].
   * Here, L is a construction parameter.
   * <p>
   * NOTE: Spikes are not allowed, all variances must be positive.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotGaussMixture : public EPScalarPotential
  {
  protected:
    // Members

    ArrayHandle<double> logp,vars; // [c_l], [v_l]
    mutable ArrayHandle<double> qbuff;
    double maxV;                   // max_l v_l
    double lseC;                   // logsumexp(c)

  public:
    // Public methods

    /**
     * Default constructor. Sets all v_l=1, all c_l=0 (meaning
     * that p_l=1/L).
     *
     * @param numl Number of mixture components (>=2)
     */
    explicit EPPotGaussMixture(int numl) {
      if (numl<2)
	throw InvalidParameterException("At least 2 components");
      logp.changeRep(numl); vars.changeRep(numl); qbuff.changeRep(numl);
      std::fill(logp.p(),logp.p()+numl,0.0);
      std::fill(vars.p(),vars.p()+numl,1.0);
      maxV=1.0;
      lseC=log((double) numl);
    }

    int numPars() const {
      return 2*vars.size();
    }

    int numConstPars() {
      return 1;
    }

    virtual double getVariance(int l) const {
      if (l<0 || l>=vars.size())
	throw OutOfRangeException(EXCEPT_MSG(""));

      return vars[l];
    }

    virtual void setVariances(const double* v) {
      int n=vars.size();
      const double* vP=std::min_element(v,v+n);

      if (*vP<(1e-16))
	throw InvalidParameterException(EXCEPT_MSG(""));
      std::copy(v,v+n,vars.p());
      vP=std::max_element(v,v+n); maxV=*vP;
    }

    virtual double getCVal(int l) const {
      if (l<0 || l>=vars.size()-1)
	throw OutOfRangeException(EXCEPT_MSG(""));

      return logp[l];
    }

    /**
     * @param cv New c_l values (size L-1)
     */
    virtual void setCVals(const double* cv) {
      int n=vars.size()-1;

      std::copy(cv,cv+n,logp.p());
      MYASS(logp[n]==0.0);
      lseC=logsumexp(logp.p(),n+1);
    }

    void getPars(double* pv) const {
      int numl=vars.size();

      pv[0]=(double) numl;
      std::copy(logp.p(),logp.p()+(numl-1),pv+1);
      std::copy(vars.p(),vars.p()+numl,pv+numl);
    }

    void setPars(const double* pv) {
      int numl=vars.size();

      if (!isValidPars(pv))
	throw InvalidParameterException(EXCEPT_MSG(""));
      setCVals(pv+1);
      setVariances(pv+numl);
    }

    /**
     * @param pv Parameter vector
     * @return   Is configuration 'pv' valid?
     */
    bool isValidPars(const double* pv) const {
      int i,numl=vars.size();

      i=(int) ceil(*pv);
      if (i!=numl || ((double) i)!=*pv)
	return false;
      const double* vP=std::min_element(pv+numl,pv+2*numl);

      return (*vP>=(1e-16));
    }

    bool suppFractional() const {
      return false;
    }

    bool isLogConcave() const {
      return false;
    }

    bool compMoments(const double* inp,double* ret,double* logz=0,
		     double eta=1.0) const;

  protected:
    // Internal methods

    /**
     * Does job of 'compMoments', but based on the unnormalized cavity
     * marginal exp( beta{-} s - 0.5 pi{-} s^2) instead of N(s|mu{-},rho{-}).
     * This means that the value log Z_hat is based on
     *   Z_hat = int t(s) exp( beta{-} s - 0.5 pi{-} s^2) d s,
     * it has to be corrected in 'compMoments'.
     * The method formally works as long as
     *   1 + pi{-} (max_l v_l) >= 1e-16.
     * In contrast, 'compMoments' requires pi{-} to be bounded away from 0.
     *
     * @param cbeta Cavity parameter beta{-}
     * @param cpi   Cavity parameter pi{-}
     * @param alpha Value alpha ret. here
     * @param nu    Value nu ret. here
     * @param logzh Value log Z_hat ret. here. Optional
     * @return      Success?
     */
    bool compMomentsInt(double cbeta,double cpi,double& alpha,double& nu,
			double* logzh) const;

    /**
     * @param a Input vector a
     * @param n Length
     * @return  log sum_k exp(a[k])
     */
    static double logsumexp(const double* a,int n);
  };

  inline bool
  EPPotGaussMixture::compMoments(const double* inp,double* ret,double* logz,
				 double eta) const
  {
    double cpi,cbeta,cmu=inp[0],crho=inp[1];
    bool rstat;

    if (eta!=1.0)
      throw NotImplemException("Fractional updates not implemented");
    if (crho<(1e-16))
      throw NumericalException(EXCEPT_MSG(""));
    cpi=1.0/crho; cbeta=cmu/crho;
    if (rstat=compMomentsInt(cbeta,cpi,ret[0],ret[1],logz)) {
      if (logz!=0)
	*logz-=0.5*(cbeta*cmu+log(crho)+SpecfunServices::m_ln2pi);
    }

    return rstat;
  }

  /*
   * Adapted from 'EPPotGaussMixture' in module 'epscal/potentials'. We use
   * z_l = 1/(1 + pi{-} v_l) here instead of rho_l
   *
   * All positive expectations are computed via logsumexp:
   *   log Z_hat
   *   log (A_til)_k  = log E_r[z_l^k]
   */
  inline bool
  EPPotGaussMixture::compMomentsInt(double cbeta,double cpi,double& alpha,
				    double& nu,double* logzh) const
  {
    int l,numl=vars.size();
    double temp,temp2,bmsq,vl,mxlz=0.0,mxla=0.0,mxla2=0.0,logz,loga,loga2;
    //bool doDebug = (fabs(cbeta)<1e-7 && fabs(cpi-1.0)<1e-7);

    // Natural parameters of "cavity" (may be undefined)
    if (1.0+cpi*maxV<(1e-16))
      return false;
    bmsq=cbeta*cbeta;
    //if (doDebug)
    //  cout << "cbeta=" << cbeta << ",cpi=" << cpi << endl;
    // First loop: log Z_l -> 'qbuff'. Maxima for accumulators
    for (l=0; l<numl; l++) {
      vl=vars[l];
      //if (doDebug)
      //cout << "l=" << l << ": v=" << vl << ", p=" << exp(logp[l]-lseC)
      //     << endl;
      temp2=-log1p(cpi*vl); // log z_l
      qbuff[l]=temp=logp[l]-lseC+0.5*(bmsq*vl/(1.0+cpi*vl)+temp2); // log Z_l
      if (l==0 || temp>mxlz) mxlz=temp;
      temp+=temp2;
      if (l==0 || temp>mxla) mxla=temp;
      temp+=temp2;
      if (l==0 || temp>mxla2) mxla2=temp;
    }
    // Second loop: Accumulation log Z_hat, log (A_til)_k, k=1,2
    logz=loga=loga2=0.0;
    for (l=0; l<numl; l++) {
      temp=qbuff[l]; // log Z_l
      temp2=-log1p(cpi*vars[l]); // log z_l
      logz+=exp(temp-mxlz);
      temp+=temp2;
      loga+=exp(temp-mxla);
      temp+=temp2;
      loga2+=exp(temp-mxla2);
    }
    // Finalize values
    logz=log(logz)+mxlz; // log Z_hat
    //if (doDebug)
    //  for (l=0; l<numl; l++)
    //cout << "r[" << l << "]=" << exp(qbuff[l]-logz) << endl;
    loga=log(loga)+mxla-logz;
    loga2=log(loga2)+mxla2-logz;
    temp=exp(loga); // A_til = E_r[z_l]
    alpha=-cbeta*temp;
    nu=temp*cpi-bmsq*exp(loga2)+alpha*alpha;
    //if (doDebug)
    //  cout << "A_til=" << temp << ", alpha=" << alpha << ", nu=" << nu
    //   << ", E[z^2]=" << exp(loga2) << endl;
    if (logzh!=0) *logzh=logz;

    return true;
  }

  /*
   * Code copied from 'FastUtils::logsumexp'
   */
  inline double EPPotGaussMixture::logsumexp(const double* a,int n)
  {
    double mx,temp,sum;

    if (n<=0) return 0.0;
    mx=*(a++); sum=1.0;
    for (int i=1; i<n; i++) {
      temp=*(a++);
      if (temp<=mx)
	sum+=exp(temp-mx);
      else {
	sum=sum*exp(mx-temp)+1.0;
	mx=temp;
      }
    }

    return mx+log(sum);
  }
//ENDNS

#endif
