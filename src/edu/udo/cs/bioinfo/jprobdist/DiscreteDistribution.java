/*
 * DiscreteDistribution.java
 *
 * Created on October 20, 2007, 10:32 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * This abstract class implements basic features of arbitrary discrete distributions.
 * The density is zero everywhere, except at the probability atoms where it is infinite.
 * The log-cdf is defined using lnP(Interval). This is an arbitrary decision.
 *
 * @author Sven Rahmann
 */
public abstract class DiscreteDistribution extends AbstractDistribution {
  
  // lnP not implemented
  // isAtom not implemented
  // isFinite not implemented
  // moment not implemented
  // cmoment not implemented
  // support not implemented
  // entropy not implemented
  // conv not implemented
   
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
    return(lnP(new Interval(min(),x)));
  }
  
  public double lnucdf(double x) {
    return(lnP(new Interval(x,max())));
  }
  
 
  /** the entropy of this discrete distribution */
  public abstract double entropy();
  // Note: If there are only finitely many atoms, one could iterate over them.
  // See FiniteDistribution.

  
  /** the log-likelihood of the given data sample under this discrete distribution */
  public double logLikelihood(double... X) {
    double ll = 0.0;
    for (double x: X) {
      ll += lnP(x);
    }
    return ll;
  }
  
}
