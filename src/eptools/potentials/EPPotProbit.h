/* -------------------------------------------------------------------
 * LHOTSE: Toolbox for adaptive statistical models
 * -------------------------------------------------------------------
 * Project source file
 * Module: eptools
 * Desc.:  Header class EPPotProbit
 * ------------------------------------------------------------------- */

#ifndef EPTOOLS_EPPOTPROBIT_H
#define EPTOOLS_EPPOTPROBIT_H

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "src/eptools/potentials/EPScalarPotential.h"
#include "src/eptools/potentials/SpecfunServices.h"
#include "src/eptools/potentials/quad/QuadPotProximalNewton.h"

//BEGINNS(eptools)
  /**
   * Probit (Gaussian c.d.f.) potential:
   *   t(s) = Phi(y (s + soff))       ['hardStep'==false]
   *   t(s) = I{y (s + soff) >= 0}    ['hardStep'==true]
   * Here, y must be in {-1,+1}.
   * Parameters: y, soff. 'hardStep' is not a parameter.
   * <p>
   * We implement 'QuadPotProximalNewton' here in order to support debugging
   * quadrature code. Implemented for 'hardStep'==false only.
   *
   * @author  Matthias Seeger
   * @version %I% %G%
   */
  class EPPotProbit : public EPScalarPotential
  {
  protected:
    // Members

    double yscal,soff;
    bool hardStep;

  public:
    // Public methods

    explicit EPPotProbit(double py,double psoff=0.0,bool phardStep=false) :
      soff(psoff),hardStep(phardStep) {
      setTarget(py);
    }

    explicit EPPotProbit(bool phardStep=false) :
      soff(0.0),hardStep(phardStep) {
      setTarget(1.0);
    }

    virtual double getTarget() const {
      return yscal;
    }

    virtual void setTarget(double py) {
      if (py!=-1.0 && py!=1.0)
	throw InvalidParameterException(EXCEPT_MSG(""));
      yscal=py;
    }

    virtual double getSOff() const {
      return soff;
    }

    virtual void setSOff(double psoff) {
      soff=psoff;
    }

    virtual bool getHardStep() const {
      return hardStep;
    }

    int numPars() const {
      return 2;
    }

    void getPars(double* pv) const {
      pv[0]=yscal; pv[1]=soff;
    }

    void setPars(const double* pv) {
      setTarget(pv[0]); setSOff(pv[1]);
    }

    bool isValidPars(const double* pv) const {
      return (pv[0]==1.0 || pv[0]==-1.0);
    }

    bool isLogConcave() const {
      return true;
    }

    bool compMoments(double cmu,double crho,double& alpha,double& nu,
		     double* logz=0,double eta=1.0) const;

    // 'QuadPotProximal' methods

    bool hasFirstDerivatives() const {
      if (hardStep) throw NotImplemException(EXCEPT_MSG(""));
      return true;
    }

    bool hasSecondDerivatives() const {
      if (hardStep) throw NotImplemException(EXCEPT_MSG(""));
      return true;
    }

    bool hasWayPoints() const {
      if (hardStep) throw NotImplemException(EXCEPT_MSG(""));
      return true;
    }

    void getInterval(double& a,bool& aInf,double& b,bool& bInf,
		     ArrayHandle<double>& wayPts) const {
      if (hardStep) throw NotImplemException(EXCEPT_MSG(""));
      aInf=bInf=true;
      wayPts.changeRep(0);
    }

    double eval(double s,double* dl=0,double* ddl=0) const {
      double temp,temp2;

      temp=yscal*(s+soff);
      temp2=-yscal*SpecfunServices::derivLogCdfNormal(temp);
      if (dl!=0) *dl=temp2;
      if (ddl!=0) *ddl=temp2*(temp2-temp*yscal);

      return -SpecfunServices::logCdfNormal(temp);
    }

    void initBracket(double h,double rho,double& l,double& r) const {
      double temp=rho*SpecfunServices::m_sqrt2/SpecfunServices::m_sqrtpi;

      l=h;
      r=yscal*std::max(0.0,yscal*(h+soff-temp))-soff;
      if (r<l) {
	temp=l; l=r; r=temp;
      }
    }
  };

  /*
   * Hyperparameter soff:
   * In log Z and all its derivatives, we replace mu{-} by mu{-} + soff, but
   * elsewhere we leave mu{-} unchanged.
   */
  inline bool
  EPPotProbit::compMoments(double cmu,double crho,double& alpha,double& nu,
			   double* logz,double eta) const
  {
    double fct,cmupbt,crhop1;

    if (eta!=1.0)
      throw NotImplemException(EXCEPT_MSG(""));
    if (crho<=0.0 || (hardStep && crho<=1e-12))
      return false;
    cmupbt=cmu+soff;
    crhop1=hardStep?crho:(crho+1.0);
    fct=yscal/sqrt(crhop1);
    alpha=cmupbt*fct;
    if (logz!=0)
      *logz=SpecfunServices::logCdfNormal(alpha); // log Z
    alpha=fct*SpecfunServices::derivLogCdfNormal(alpha); // alpha
    nu=alpha*(alpha+cmupbt/crhop1); // nu
    // DEBUG!!
    //if (fabs(cmu-(-11.276141))<1e-5)
    //  cout << "DEBUG: y=" << yscal << ",soff=" << soff << "hStep="
    //   << hardStep << ",cmu=" << cmu << ",crho=" << crho
    //   << ",alpha=" << alpha << ",nu=" << nu << endl;
    // END DEBUG

    return true;
  }
//ENDNS

#endif
