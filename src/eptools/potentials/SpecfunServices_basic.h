class SpecfunServices
{
public:
  // Constants

  static const double M_LN2PI = 1.837877066409345339081937709125;

  // Static methods

  /**
   * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
   * log Phi(z).
   *
   * @param z Argument
   * @return  log Phi(z)
   */
  static double logCdfNormal(double z) {
    // TODO: Implement!
    throw NotImplemException(EXCEPT_MSG(""));
  }

  /**
   * If Phi(z) denotes the c.d.f. of N(0,1), this method computes
   *   f(z) = (d/dz) log Phi(z) = N(z)/Phi(z).
   * The function f(-z) is called hazard function in statistics.
   *
   * @param z Argument
   * @return  (d/dz) log Phi(z)
   */
  static double derivLogCdfNormal(double z) {
    // TODO: Implement!
    throw NotImplemException(EXCEPT_MSG(""));
  }

  /**
   * Computes natural log of Gamma(z) for z>0. Note that if z is a
   * natural number, then z! = Gamma(z+1).
   *
   * @param z Argument (positive)
   * @return  log Gamma(z)
   */
  static double logGamma(double z) {
    // TODO: Implement!
    throw NotImplemException(EXCEPT_MSG(""));
  }
};
