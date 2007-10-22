/*
 * Interval.java
 *
 * Created on October 20, 2007, 6:37 PM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;
import static java.lang.Math.*;

/**
 * This class represents an unmodifiable one-dimensional interval.
 * @author Sven Rahmann
 */
public class Interval {
  
  /** left boundary */
  public final double a;
  /** right boundary */
  public final double b;
  /** type of interval (Closed, Open, OpenClosed, ClosedOpen, NaN) */
  public final Type type;
  /** true iff this interval represents a single point */
  public final boolean isPoint;
  /** true iff this interval is empty (and not NaN) */
  public final boolean isEmpty;
 
  
  /** create a new interval (a,b), [a,b], (a,b], or [a,b),
   * depending on the given Type t
   *@param a  left boundary
   *@param b  right boundary
   *@param type  interval type (e.g., Inteval.Type.Closed)
   */
  public Interval(double a, double b, Type type) {
    this.a=a;
    this.b=b;
    this.type=(Double.isNaN(a) || Double.isNaN(b))? Type.NaN : type;
    this.isPoint = (b==a && type==Type.Closed);
    this.isEmpty = (b<a || (b==a && type!=Type.Closed));
  }
  
  /** create a new closed interval [a,b] */
  public Interval(double a, double b) {
    this(a,b, Type.Closed);
  }
  
  /** create a new point inteval [b,b] */
  public Interval(double b) {
    this(b, b, Type.Closed);
  }
 
  /** types of one-dimensional intervals */
  public enum Type {
    /** represents an interval (a,b) that is open at both boundaries */
    Open, 
    /** represents an interval [a,b] that is closed at both boundaries */
    Closed, 
    /** represents an interval (a,b] that is open at the left and closed at the right boundary */
    OpenClosed, 
    /** represents an interval [a,b) that is closed at the left and open at the right boundary */
    ClosedOpen,
    /** represents an interval that is not well-defined (one of its boundaries is NaN) */
    NaN,
  }
  
  /** return a closed(!) interval whose endpoints are integer-divisible by eps
   * and that is completely contained in the present interval.
   *@param eps  the granularity (e.g., 1.0 for integer intervals)
   *@return contained closed interval, or a NaN interval
   */
  public Interval getContainedEpsInterval(final double eps) {
    if (eps==0) return this;
    double f = eps*ceil(this.a/eps);
    double l = eps*floor(this.b/eps);
    switch(this.type) {
      case Open:
        if (f==this.a) f+=eps;
      case ClosedOpen:
        if (l==this.b) f-=eps;
      case Closed:
        break;
      case OpenClosed:
        if (f==this.a) f+=eps;
        break;
      default:
        f = Double.NaN;
        l = Double.NaN;
        break;
    }
    return new Interval(f,l, (this.type!=type.NaN? Type.Closed: Type.NaN));
  }

  
  // Finally, we define some useful constants:
  
  /** the empty interval */
  public static final Interval AnEmptyInterval = new Interval(0,0,Type.Open);
  /** the whole real line (open at infinity) */
  public static final Interval TheRealLine      = new Interval(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, Type.Open);
  
}

