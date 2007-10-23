/*
 * BinomialDistribution.java
 *
 * Created on October 21, 2007, 1:33 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;
import java.util.Iterator;
import static java.lang.Math.*;
import static edu.udo.cs.bioinfo.jprobdist.MathFunctions.*;
import java.util.NoSuchElementException;

/**
 * This class implements the Binomial distribution, i.e.,
 * the distribution of the number of successes in n independent trials,
 * where each trial has a success probability of p.
 *
 * @author Sven Rahmann
 */
public class BinomialDistribution extends FiniteDistribution
    implements UnimodalDistribution {
  
  /** the number of trials */
  public final long n;
  /** the natural logarithm of the success probability */
  public final double lnp;
  /** the natural logarithm of the failure probability */
  public final double lnq;
  /** the expectation of this distribution */
  public final double E;
  
  /**
   * Generates a new Binomial Distribution with parameters n and p.
   * @param n The nonnegative number of trials (has to be >= 0)
   * @param p The probability for a success
   */
  public BinomialDistribution(final long n, final double p) {
    if (n<0) throw new IllegalArgumentException("n must be nonnegative");
    if (p<0 || p>1) throw new IllegalArgumentException("p must be in [0,1]");
    this.n = n;
    this.lnp = log(p);
    this.lnq = log1p(-p);
    this.E = n*p;
  }
  
  /**
   * Generates a new Binomial Distribution with parameters n and p,
   * where p is approximately equal to 1, and one wants to give
   * the failure parameter q explicitly.
   * @param n  The nonnegative number of trials (has to be >= 0)
   * @param p  The probability for a success, must be given as 1.0 in this constructor,
   * but is calculated as 1-q.
   *@param q The probability for a failure (must be between 0 and 1).
   *
   */
  public BinomialDistribution(final long n, final double p, final double q) {
    if (n<0) throw new IllegalArgumentException("n must be nonnegative");
    if (p!=1.0) throw new IllegalArgumentException("p must be 1.0 if you want to specify q");
    if (q<0 || q>1) throw new IllegalArgumentException("q must be in [0,1]");
    this.n = n;
    this.lnp = log1p(-q);
    this.lnq = log(q);
    this.E = exp(log(n)+lnp);
  }
  

  public double lnP(double x) {
    if (!isAtom(x)) return Double.NEGATIVE_INFINITY;
    if (x==0) return ( (lnp==Double.NEGATIVE_INFINITY)? 0.0 : n*lnq);
    if (x==n) return ( (lnq==Double.NEGATIVE_INFINITY)? 0.0 : n*lnp);
    return(lnbincoeff(n,x)+x*lnp+(n-x)*lnq);
  }
  
  
  @Override
  public double E() {
    return E;
  }
  
  @Override
  public double Var() {
    return E*exp(lnq);
  }
  
  @Override
  public double skewness() {
    return (1-2*exp(lnp))/std();
  }
  
  @Override
  public double kurtosisExcess() {
    return (1-6*exp(lnp+lnq))/Var();
  }
  
  
  public final double min() {
    return 0;
  }
  
  public final double max() {
    return n;
  }
  
  public double median() {
    // the median is either floor(E), or one less, or one more [?]
    final double x = floor(E)-1;
    assert(cdf(x-1)<0.5);
    if (cdf(x)>=0.5) return x;
    if (cdf(x+1)>=0.5) return x+1;
    if (cdf(x+2)>=0.5) return x+2;
    throw new AssertionError(String.format("Unexpected: cdf(x+2)=%f < 0.5",cdf(x+2)));
  }
  
  public final double mode() {
    return ceil(exp(lnp)*(n+1))-1;
    // this is the smallest mode
  }
  
  public final Interval modeInterval() {
    return new Interval(mode(), floor(exp(lnp)*n), Interval.Type.Closed);
  }
  
  
  public boolean isAtom(final double x) {
    return (x>=0 && x<=n && x==MathFunctions.xround(x));
  }
  
  public final double closestAtom(final double x) {
    if (x<=0) return 0;
    if (x>=n) return n;
    return xround(x);
  }
  
  public Iterator<Double> iterator() {
    return new LatticeIterator(n, 0.0, 1.0);
  }
  
  public Iterator<Double> iterator(Interval ab) {
    final Interval fl = ab.getContainedEpsInterval(1.0);
    final long first = (fl.a<=0)? 0 : (long)fl.a;
    final long last  = (fl.b>=n)? n : (long)fl.b;
    return new LatticeIterator(last-first, first, 1.0);
  }

  
  /** Estimate the Binomial success parameter p from a data sample X,
   * whose elements are assumed to be realization of Binomial(n,p)
   * with known n and unknown p.
   * The ML estimator is given by sum(X)/(n*X.length).
   *@param n the assumed numer of trials for each element of X
   *@param X the observed number of successes
   *@return the ML estimate of the success probability
   */
  public static final double pFromSample(final long n, final double... X) {
    double sum = 0.0;
    for (double x: X) sum+=x;
    return sum / (n*X.length);
  }
  
}
