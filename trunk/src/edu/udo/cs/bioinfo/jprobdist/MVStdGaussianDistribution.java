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

  public double f(double... x) {
    if (x.length!=3) throw new IllegalArgumentException();
    return 1.0;
  }

  public double P(double... x) {
    return 1.0;
  }


  public double cdf(final double... x) {
    return 1.0;
  }

  public double ucdf(final double... x) {
    return 1.0;
  }

  public boolean isAtom(final double... x) {
    return false;
  }
  
  public boolean isFinite() {
    return false;
  }
  
  
}
