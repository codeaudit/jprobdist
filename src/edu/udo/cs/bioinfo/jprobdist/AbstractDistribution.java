/*
 * AbstractDistribution.java
 *
 * Created on October 20, 2007, 9:18 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * This class represents an abstract implementation of a
 * univariate (i.e., 1-dimensional) probability distribution.
 *
 * It contains default implementations of frequently used convenience
 * methods.
 *
 * For concrete probability distribution, more efficient methods
 * may exist, so the methods in this class should be overrided where desired.
 *
 * @author Sven Rahmann
 */
public abstract class AbstractDistribution implements UVDistribution {
  
  // lnf not implemented
  // lnP not implemented
  // lncdf not implemented
  // lnucdf not implemented
  // isAtom not implemented
  // isFinite not implemented
  // moment, cmoment not implemented
  // support, min, max not implemented
  // conv not implemented
  
  /* NOTE: P(double x), P(Interval ab)
   * return the probability mass at point x resp. in the interval ab.
   *
   * In principle, there are three ways to compute this probability,
   * depending on the location of the interval.
   * 1. We can base the computation on the cdf,
   * if no cancellation problems occur when computing cdf(b) - cdf(a),
   * and if cdf() has an efficient implementation.
   * 2. We can base the computation on the upper cdf,
   * if no cancellation problems occur when computing ucdf(a) - ucdf(b),
   * and if ucdf() has an efficient implementation (NOT by default!)
   * 3. We can base the computation on the sum of atomic probabilities P(x),
   * where x is in the interval ab.
   * This is preferable when the distribution is discrete and no efficient
   * cdf is given.
  public double P(Interval ab) {
    if (ab.type == Interval.Type.NaN) return Double.NaN;
    if (ab.isEmpty) return 0.0;
    if (ab.isPoint) return P(ab.a);
    switch(ab.type) { // no breaks are required because all braches return.
      case Closed:
        return(cdf(ab.b)-cdf(ab.a)+P(ab.a));
      case Open:
        return(cdf(ab.b)-cdf(ab.a)-P(ab.b));
      case OpenClosed:
        return(cdf(ab.b)-cdf(ab.a));
      case ClosedOpen:
        return(cdf(ab.b)-cdf(ab.a)+P(ab.a)-P(ab.b));
      default:
        return Double.NaN;
    }
  }
   */
  
  // =================================================================
  // provide a default implementation of f, P, cdf, ucdf
  // by taking the exponentials of their logarithmic implementations.
  
  public double f(final double x) {
    return Math.exp(lnf(x));
  }
  
  public double P(final double x) {
    return Math.exp(lnP(x));
  }
  
  public double P(final Interval ab) {
    return Math.exp(lnP(ab));
  }
  
  public double cdf(final double x) {
    return Math.exp(lncdf(x));
  }
  
  public double ucdf(final double x) {
    return Math.exp(lnucdf(x));
  }
  
  // =====================================================================
  // provide generic implementations of E, Var, std, kurtosis, median, iqr,
  // and boxPlotStatistics
  
  public double E() {
    return moment(1.0);
  }
  
  public double Var() {
    return cmoment(2.0);
  }
  
  public double std() {
    return Math.sqrt(Var());
  }
  
  public double skewness() {
    return cmoment(3.0)/Math.pow(Var(),1.5);
  }
  
  public double kurtosisExcess() {
    final double v = Var();
    return cmoment(4.0)/(v*v) - 3.0;  // -3 so the Gaussian has zero kurtosis
  }
  
  public double kurtosisProper() {
    return kurtosisExcess() + 3.0;
  }
  
  public double median() {
    return(qf(0.5));
  }
  
  public double iqr() {
    return(qf(0.75)-qf(0.25));
  }
  
  
  public final double[] boxPlotStatistics() {
    double med = median();
    double lq  = qf(0.25);
    double uq  = qf(0.75);
    double lw  = lq - 1.5*iqr();
    double uw  = uq + 1.5*iqr();
    return new double[] {lw, lq, med, uw, uq};
  }
  
  // ======================================================================
  // provide a generic quantile function
  // implemented by a simple bisection method
  // qf(p) := inf {x : cdf(x) >= p}
  
  public double qf(final double p) {
    if (p<0 || p>1)
      throw new IllegalArgumentException("qf(p): p must be in [0,1], is "+Double.toString(p));
    if (p==0.0) return Double.NEGATIVE_INFINITY;
    if (p==1.0) return max();
    // now 0<p<1
    
    // find an initial interval [L,R] such that p is in [cdf(L),cdf(R)]
    double L = min();
    if (L==Double.NEGATIVE_INFINITY) {
      L=Double.MIN_VALUE;
    }
    double R = max();
    if (R==Double.POSITIVE_INFINITY) {
      R=Double.MAX_VALUE;
    }
    
    // loop until L==R
    while (L!=R) {
      final double C    = (L<0 && R>0)? (R/2-L/2)+L : (L/2+R/2);
      final double Ccdf = cdf(C);
      if (p<=Ccdf) R=C; else L=C;
    }
    
    assert(cdf(L)>=p); // if this ever fails, we need to write more careful code
    return L;
  }
  
  
  // =======================================================================
  // provide default random() method by using the inversion method,
  // i.e. calling qf on a uniform random variable, from Math.random()
  
  public double random() {
    final double p = Math.random();
    if (p==0.0) return max();
    return qf(p);
  }
  
  public double[] random(final int n) {
    final double[] r = new double[n];
    for(int i=0; i<n; i++) r[i]=random();
    return r;
  }
  
  
  // ======================================================
  // Delegate all vararg functions to univariate functions
  // ======================================================
  
  public final double f(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return f(x[0]);
  }
  public final double lnf(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return lnf(x[0]);
  }
  
  public final double P(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return P(x[0]);
  }
  public final double lnP(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return lnP(x[0]);
  }
  
  public final double cdf(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return cdf(x[0]);
  }
  public final double lncdf(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return lncdf(x[0]);
  }
  
  public final double ucdf(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return ucdf(x[0]);
  }
  public final double lnucdf(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return lnucdf(x[0]);
  }
  
  public final boolean isAtom(final double... x) {
    if (x.length!=1) throw new DimensionMismatchException();
    return isAtom(x[0]);
  }
  
  
}
