/*
 * DiscreteDistribution.java
 *
 * Created on October 20, 2007, 10:32 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

import java.util.Iterator;
import java.util.NoSuchElementException;
import static java.lang.Math.*;

/**
 * This abstract class implements basic features of arbitrary discrete distributions.
 * The density is zero everywhere, except at the probability atoms where it is infinite.
 * The log-cdf is defined using lnP(Interval). This is an arbitrary decision.
 *
 * @author Sven Rahmann
 */
public abstract class DiscreteDistribution extends AbstractDistribution
    implements Iterable<Double> {
  
  // abstract:
  // lnP
  // isAtom
  // isFinite
  // E(h())
  
  public final double f(double x) {
    // For discrete distributions, the density is zero, except at atoms,
    // where it's infinite
    return (isAtom(x)? Double.POSITIVE_INFINITY : 0.0);
  }
  
  public final double lnf(double x) {
    // For discrete distributions, the log-density is -inf, except at atoms,
    // where it's +inf
    return (isAtom(x)? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY);
  }
  
  public double lncdf(double x) {
    // For discrete distributions, the log-cdf might be conveniently computed from lnP.
    // But this should be overridden for efficiency, where possible.
    if(x>=max()) return 0.0;
    return(lnP(new Interval(min(),x)));
  }
  
  public double lnucdf(double x) {
    if(x>=min()) return 0.0;
    return(lnP(new Interval(x,max())));
  }
  
  // Since the distribution is discrete, we can round the result returned by qf
  // towards the nearest atom, unless it's -inf.
  public double qf(double x) {
    double result = super.qf(x);
    return (result==Double.NEGATIVE_INFINITY)? result : closestAtom(result);
  }
  
  /** the expectation of a real-valued function under this probability distribution */
  public double E(MathFunctions.RealFunction h) {
    double r = 0.0;
    double ptotal = 0.0;
    int pdec = 0;
    int rsame = 0;
    double p=0, op=0;
    double term = Double.NEGATIVE_INFINITY, oterm = term; // magnitude of old term
    
    for(double x : this) { // iterate over all atoms. Could be infinitely many.
      op = p; oterm = term;
      ptotal += (p = P(x));
      if(p==0.0) continue;
      if (p<op) pdec++; else pdec=0;
      final double hx = h.valueAt(x);
      term = p*hx;
      if (r+term==r) rsame++; else rsame=0;
      r += term; 
      if ( (abs(p-1.0)<16*MathFunctions.DBL_TOL) && (pdec>=20 && rsame>=20) ) break;
    }
    return r;
  }
 
  
  /** The expectation itself */
  public double E() {
    return E(MathFunctions.ID);
  }
  
  /** The entropy of this discrete distribution.
   *  It is defined as E[-ln(P)] under P.
   */
  public double entropy() {
    final MathFunctions.RealFunction h = new MathFunctions.RealFunction() {
      public double valueAt(double x) { return -lnP(x); }
    };
    return E(h);
  }
  
  /** The m-th moment of this distribution.
   *  It is defined as E[X^m].
   */
  public double moment(final double m) {
    final MathFunctions.RealFunction h = new MathFunctions.RealFunction() {
      public double valueAt(double x) { return Math.pow(x,m); }
    };
    return E(h);
  }
  
  /** The m-th central moment of this distribution.
   *  It is defined as E[(X-EX)^m]
   */
  public double cmoment(final double m) {
    final double ee = E();
    final MathFunctions.RealFunction h = new MathFunctions.RealFunction() {
      public double valueAt(double x) { return Math.pow(x-ee,m); }
    };
    return E(h);
  }
  
 
  
  /** the log-likelihood of the given data sample under this discrete distribution */
  public double logLikelihood(double... X) {
    double ll = 0.0;
    for (double x: X) {
      ll += lnP(x);
    }
    return ll;
  }
  
  
  /** returns an iterator that iterates over all atoms of the
   * present finite distribution
   */
  public abstract Iterator<Double> iterator();
  
  /** returns an iterator that iterates over the atoms of
   *the present finite distribtion in the given interval
   */
  public abstract Iterator<Double> iterator(Interval ab);
  
  
  
  /** a general purpose discrete iterator:
   * It iterates through the n+1 values (shift + scale*[0..n]).
   *  Using n=Integer.MAX_VALUE will result in an infinite loop.
   */
  public static class LatticeIterator implements Iterator<Double> {
    private long i;
    private final long n;
    private final double lambda;
    private final double mu;
    
    /** iterates through the n+1 values shift + scale*[0..n].
     *  Using n=Integer.MAX_VALUE will result in an infinite loop.
     */
    public LatticeIterator(final long n, double shift, double scale) {
      i=0;
      this.n = n;
      lambda = scale;
      mu = shift;
    }
    
    public boolean hasNext() {
      return (i<=n);
    }
    
    public Double next() {
      if(!hasNext()) throw new NoSuchElementException(String.valueOf(i));
      return(Double.valueOf(mu+(lambda*i++)));
    }
    
    public void remove() {
      throw new UnsupportedOperationException();
    }
  } // end of iterator class
  
  
}
