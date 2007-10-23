/*
 * ContinuousDistribution.java
 *
 * Created on October 21, 2007, 1:09 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * This abstract class implements basic features of arbitrary continuous distributions.
 * The probability mass is zero everywhere, and there are no atoms. 

 * @author Sven Rahmann
 */
public abstract class ContinuousDistribution extends AbstractDistribution {
 
  // f not implemented
  // P(Interval) not implemented
  // cdf, ucdf not implemented
  // moment not implemented
  // cmoment not implemented
  // support not implemented
  // conv not implemented
  

  public final double P(final double x) {
    return 0.0;
  }
  
  public final boolean isAtom(final double x) {
    return false;
  }
  
  public final double closestAtom(final double x) {
    return Double.NaN;
  }
  
  public final boolean isFinite() {
    return false;
  }
 
  
}
