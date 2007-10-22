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
public abstract class FiniteDistribution extends DiscreteDistribution 
    implements Iterable<Double> {

  // does not implement lnP
  // does not implement min
  // does not implement max
  // does not implement iterator


  // explicit moment implementation by iterating over all atoms
  public double moment(final double m) {
    double r = 0.0;
    for(double x : this) {
      final double lp = lnP(x);
      r += Math.exp(lp) * Math.pow(x, m);
    }
    return r;
  }

  // explicit central moment implementation by iterating over all atoms
  public double cmoment(final double c) {
    final double E = E();
    double r = 0.0;
    for(double x : this) {
      final double lp = lnP(x);
      r += Math.exp(lp) * Math.pow(x-E, c);
    }
    return r;
  }

  
  public Interval support() {
    return new Interval(min(),max(),Interval.Type.Closed);
  }

  public final boolean isFinite() {
    return true;
  }

  // an explicit implementation of the entropy function by summation
  public double entropy() {
    double H = 0.0;
    for(double x : this) {
      final double lp = lnP(x);
      H += lp * Math.exp(lp);
    }
    assert(H<=0);
    return -H;
  }
  
  // an inefficient implementation of testing whether x is an atom
  public boolean isAtom(double x) {
    for(double xx : this) {
      if(xx==x) return true;
    }
    return false;
  }
  
  /** returns an iterator that iterates over all atoms of the
   * present finite distribution
   */
  public abstract Iterator<Double> iterator();
  
  /** returns an iterator that iterates over the atoms of
   *the present finite distribtion in the given interval
   */
  public abstract Iterator<Double> iterator(Interval ab);
  
}
