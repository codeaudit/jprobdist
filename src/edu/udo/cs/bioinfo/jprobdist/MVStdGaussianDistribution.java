/*
 * MVStdGaussianDistribution.java
 *
 * Created on October 20, 2007, 5:47 PM
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * a multivariate Standard Gaussian Distribution
 * For DEMONSTRATION PURPOSES ONLY
 *
 * @author Sven Rahmann
 */
public class MVStdGaussianDistribution implements Distribution {
  
  /** number of dimensions */
  private final int d; 
  
  /** Creates a new instance of MVStdGaussianDistribution */
  /*@param d  number of dimensions
   */
  public MVStdGaussianDistribution(int d) {
    this.d = d;
  }

  public double lnf(double... x) {
    return 0.0;
  }
  

  public double lnP(double... x) {
    return Double.NEGATIVE_INFINITY;
  }


  public double lncdf(final double... x) {
    return 0.0;
  }

  public double lnucdf(final double... x) {
    return 0.0;
  }

  public boolean isAtom(final double... x) {
    return false;
  }
  
  public boolean isFinite() {
    return false;
  }

 
  public double cdf(final double... x) {
    return 1.0;
  }

  public double ucdf(final double... x) {
    return 0.0;
  }

  public double f(final double... x) {
    return 1.0;
  }

  public double P(final double... x) {
    return 0.0;
  }
  
  
}
