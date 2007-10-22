/*
 * Distribution.java
 *
 * Created on October 20, 2007, 5:36 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * interface for general (i.e., univariate, bivariate or multivariate
 * discrete or continuous or mixed) distributions.
 * While the interface functionality primarily aims at probability distributions
 * (i.e., unsigned measures with total mass 1), this is never enforced.
 *
 * <p>
 * This code is intended for teaching and research purposes only.
 * This code is licensed under the Artistic License 2.0.
 * </p>
 * 
 * In the general multivariate case, points at which funcitons should be
 * evaluated are given as double arrays, or a list of doubles, using
 * Java 5's varargs mechanism.
 *
 * @author Sven Rahmann
 */
public interface Distribution {

  // ******************** Descriptors ******************************/
  
  /** density function 
   *@param x  the point x at which to compute the density f
   *@return density value f(x)
   */
  public double f(final double... x);
  
  /** logarithmized density function.
   * This function may be more efficient or more accurate than
   * calling Math.log(f(x)) explicitly.
   *@param x the point x at which to compute the log-density
   *@return ln(f(x))
   */
  public double lnf(final double... x);
  
  /** probability mass function at a single point
   *@param x the point at which to compute the probability mass
   *@return P(x), the probability mass at x
   */
  public double P(final double... x);
  
  /** logarithmized probability mass funciton at a single point 
   * This function may be more efficient or more accurate than
   * calling Math.log(P(x)) explicitly.
   *@param x the point at which to compute the log-probability mass
   *@return ln(P(x))
   */
  public double lnP(final double... x);
  
  /** cumulative distribution function
   *@param x
   *@return total mass at points <= x
   */
  public double cdf(final double... x);
  
  /** logarithm of cdf at x 
   * This function may be more efficient or more accurate than
   * calling Math.log(cdf(x)) explicitly.
   *@param x
   *@return ln(cdf(x))
   */
  public double lncdf(final double... x);
  
  /** upper cumulative distribution function 
   *@return total mass at points >= x (in all coordinates)
   */
  public double ucdf(final double... x);
  
  /** logarithm of ucdf (upper cumulative distribution function) at x,
   * where ucdf(x) = Prob(X>=x in all coordinates)
   * This function may be more efficient or more accurate than
   * calling Math.log(ucdf(x)) explicitly.
   *@param x
   *@return ln(ucdf(x))
   */
  public double lnucdf(final double... x); 
  
  /** returns true if and only if point x is an atom */
  public boolean isAtom(final double... x);
  
  /** returns true iff this is a finite distribution */
  public boolean isFinite();

  
  // ************************** Actions ******************************/
  
  // none defined
  
}
