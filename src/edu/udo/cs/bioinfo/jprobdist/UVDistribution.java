/*
 * UVDistribution.java
 *
 * Created on October 20, 2007, 6:25 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * interface for a general univariate (i.e., 1-dimensional distribution)
 * @author Sven Rahmann
 */
public interface UVDistribution extends Distribution {
  
  
  // Specialize the functions in Distribution to a 'double' argument.
  // The varargs (double... x) functions will be delegated to these functions
  // in the abstract class AbstractDistribution, which provides some 
  // default implementations for univariate distributions.
  
  public double f(final double x);  
  public double lnf(final double x);
  
  public double P(final double x);
  public double lnP(final double x);
      
  public double cdf(final double x);
  public double lncdf(final double x);
  
  public double ucdf(final double x);
  public double lnucdf(final double x);
   
  public boolean isAtom(final double x);

  
  // Now, define new functions that are new to UVDistribution

  /** probability mass in the interval ab
   *@param ab  the interval
   *@return total probability mass in ab
   */
  public double P(final Interval ab);    // new for UVDistribution
  
  /** logarithm of probability mass in the interval ab;
   * in certain cases, this function may be more efficient or more accurate
   * than calling Math.log(P(ab)).
   *@param ab  the interval
   *@return log(P(ab))
   */
  public double lnP(final Interval ab);  // new for UVDistribution

  /** expectation */
  public double E();
  
  /** variance */
  public double Var();
  
  /** standard deviation */
  public double std();
  
  /** skewness */
  public double skewness();
  
  /** The kurtosis excess of the distribution.
   * The standard Gaussian distribution has kurotsis excess zero.
   */
  public double kurtosisExcess();
  
  /** The proper kurtosis of the distribution.
   * The standard Gaussian distribution has a propoer kurtosis of +3.
   * Usually, one considers the kurtosis excess (see <code>kurtosisExcess</code>),
   * where one subtracts 3 from the proper kurtosis.
   */
  public double kurtosisProper();

  /** m-th moment 
   * @param m moment to compute 
   */
  public double moment(final double m);
  
  /** m-th central moment
   * @param c central moment to compute
   */
  public double cmoment(final double c);
  
  /** median */
  public double median();
  
  /** inter-quartile range */
  public double iqr();
  
  /** maximum (supremum) value with nonzero probability or density */
  public double max();
  
  /** minimum (infimum) value with nonzero probability or density */
  public double min();
  
  /** the interval where this probability distribution exists */
  public Interval support();
  
  /** quantile function (the inverse of the cdf)
   *@param p  a probability in [0,1]
   *@return the p-quantile, i.e. the value of x where cdf(x)=p.
   * In non-continuous distributions, the relationship between cdf and qf is more complicated. 
   * It is likely that there is no such value x that the CDF of x gives the desired probability.
   * Several possiblities of defining the qf actually exist.
   * qf(p) in that case returns the smallest x such that cdf(x)>=p, more precisely:
   * qf(p) := inf { x : cdf(x)>= p }
   *
   */
  public double qf(final double p);
  
  /** markers for Box-Whisker plot  */
  public double[] boxPlotStatistics();
  

  // ************************** Actions ******************************/

  /** convolution with another Distribution */
  //public UVDistribution conv(UVDistribution other);
  
  /** generate a random number from this distribution */
  public double   random();
  
  /** generate n random numbers from this distribution */
  public double[] random(final int n);
    
  
}
