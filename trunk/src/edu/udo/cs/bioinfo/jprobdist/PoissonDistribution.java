/*
 * PoissonDistribution.java
 *
 * Created on October 22, 2007, 8:31 PM
 */

package edu.udo.cs.bioinfo.jprobdist;

import java.util.Iterator;
import static java.lang.Math.*;
import edu.udo.cs.bioinfo.jprobdist.DiscreteDistribution.LatticeIterator;


/**
 * This class implements the Poisson Distribution
 * @author Sven Rahmann
 */
public class PoissonDistribution extends DiscreteDistribution
    implements UnimodalDistribution {
  
  /** the Poisson parameter */
  public final double lambda;   // since it's final, it can be public!
  /** the logarithm of the Poisson parameter */
  public final double lnlambda;
  ///** exp(-lambda) */
  //public final double expminuslambda;
  
  /** create a new Poisson distribution with the given parameter lambda.
   * The Poisson distribution models rare events.
   * Both expectation and variance are equal to lambda.
   */
  public PoissonDistribution(final double lambda) {
    if(lambda<0 || Double.isNaN(lambda))
      throw new IllegalArgumentException("Poisson requires lambda>0");
    this.lambda = lambda;
    lnlambda = Math.log(lambda);
    //expminuslambda = Math.exp(-lambda);
  }
  
  //==============================================================
  // implement lnP(x). P(x) is covered by superclasses.
  // Also implement P and lnP for intervals.
  
  public final double lnP(final double x) {
    if(!isAtom(x)) return Double.NEGATIVE_INFINITY;
    return(-lambda + x*lnlambda - MathFunctions.lnfactorial(x));
  }
  
  
  public final double P(final Interval ab) {
    if(ab.isEmpty) return Double.NEGATIVE_INFINITY;
    if(ab.isPoint) return lnP(ab.a);
    final double a = ab.a;
    final double b = ab.b;
    if(b-a<=10) {
      double sum=0;
      for(Iterator<Double> it=iterator(ab); it.hasNext(); sum+=it.next());
      return sum;
    }
    if (a > lambda) return (ucdf(a)-ucdf(b)+P(b));
    return (cdf(b)-cdf(a)+P(a));
  }
  
  public final double lnP(final Interval ab) {
    return log(P(ab));
  }
  
  
  // =============================================================
  // reimplement cdf, lncdf, ucdf, lnucdf
  
  public final double cdf(final double x) {
    if(lambda==0) return (x>=0? 1.0 : 0.0);
    return MathFunctions.gammaQ(Math.floor(x+1), lambda);
  }
  
  public final double lncdf(final double x) {
    return(log(cdf(x)));
  }
  
  public final double ucdf(final double x) {
    if(lambda==0) return (x<=0? 1.0 : 0.0);
    return MathFunctions.gammaP(Math.ceil(x+1), lambda);
    //TODO: Check!
  }
  
  public final double lnucdf(final double x) {
    return(log(ucdf(x)));
  }
  
  
  // =============================================================
  // implement special moments.
  
  public final double E() {
    return lambda;
  }
  
  public final double Var() {
    return lambda;
  }
  
  public final double skewness() {
    return 1/Math.sqrt(lambda);
  }
  
  public final double kurtosisExcess() {
    return 1/lambda;
  }
  
  // ============================================================
  // implement simple stuff
  public final double max() {
    return (lambda==0)? 0 : Double.POSITIVE_INFINITY;
  }
  
  public final double min() {
    return 0.0;
  }
  
  public final Interval support() {
    return (lambda==0)? Interval.TheZeroInterval : Interval.TheNonnegativeReals;
  }
  
  public final boolean isFinite() {
    return (lambda==0);
  }
  
  public boolean isAtom(final double x) {
    return (x>=0 && x==MathFunctions.xround(x));
  }
  
  public final double closestAtom(final double x) {
    if(x<=0) return 0;
    return (MathFunctions.xround(x));
  }
  
  // ================================================================
  // implement mode functions
  // the mode is floor(lambda). If lambda is an integer, lambda-1 is also a mode.
  // mode() is supposed to return the smallest mode.
  public final double mode() {
    if (lambda==0) return 0;
    if(MathFunctions.xround(lambda)==lambda) return lambda-1;
    return Math.floor(lambda);
  }
  
  public final Interval modeInterval() {
    if (lambda==0) return Interval.TheZeroInterval;
    if(MathFunctions.xround(lambda)==lambda) return new Interval(lambda-1,lambda);
    return new Interval(lambda);
  }
  
  
  // ===================================================================
  // implement atom iterators: 0, 1, 2, ... ad "infinitum"
  public Iterator<Double> iterator() {
    return new DiscreteDistribution.LatticeIterator(Long.MAX_VALUE, 0.0, 1.0);
  }
  
  // iterate over integers in interval ab.
  public Iterator<Double> iterator(Interval ab) {
    final Interval fl = ab.getContainedEpsInterval(1.0);
    final long first = (fl.a<=0)? 0 : (long)fl.a;
    final long last  = (long)fl.b;
    return new LatticeIterator(last-first, first, 1.0);
  }
  
  
}
