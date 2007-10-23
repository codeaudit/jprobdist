/*
 * FiniteDistribution.java
 *
 * Created on October 21, 2007, 5:03 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;
import java.util.Iterator;

/**
 * This class provides basic support for arbitrary finite probability distributions.
 * @author Sven Rahmann
 */
public abstract class FiniteDistribution extends DiscreteDistribution {

  // Missing implementations for finite distributions:
  // lnP(x)
  // min(), max()
  // iterator

  
  // simply add up the finitely many atom probabilities in the interval.
  // Should be overridden with more efficient methods in subclasses.
  public double lnP(final Interval ab) {
    double p = Double.NEGATIVE_INFINITY;
    Iterator<Double> it = iterator(ab);
    while (it.hasNext()) p = MathFunctions.logsum(p, lnP(it.next()));
    return p;
  }
  
  // simply add up the finitely many atom probabilities in the interval.
  // Should be overridden with more efficient methods in subclasses.
  @Override
  public double P(final Interval ab) {
    double p = 0.0;
    Iterator<Double> it = iterator(ab);
    while (it.hasNext()) p += P(it.next());
    return p;
  }  
  
  // enumerating implementation of arbitrary Expectation function.
  // Cannot be overridden in subclasses.
  // But particular expectation functions,
  // such as E(), moment(m), cmoment(m), entropy(),
  // should be overridden where possible for efficiency!
  public final double E(final MathFunctions.RealFunction h) {
    double r = 0.0;
    for(double x : this) {
      final double p = P(x);
      if (p==0.0) continue;
      final double hx = h.valueAt(x);
      if (hx==0.0) continue;
      r += p * hx;
    }
    return r;
  }

   
  // an inefficient implementation of testing whether x is an atom.
  // Should be overriden in subclasses!
  public boolean isAtom(double x) {
    for(double xx : this) {
      if(xx==x) return true;
    }
    return false;
  }

  // the support of a finite distribution is always closed and compact.
  public Interval support() {
    return new Interval(min(),max(),Interval.Type.Closed);
  }

  // of course it's finite. Cannot be overriden in subclasses.
  public final boolean isFinite() {
    return true;
  }

  
}
