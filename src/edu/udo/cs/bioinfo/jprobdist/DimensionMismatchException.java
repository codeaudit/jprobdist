/*
 * DimensionMismatchException.java
 *
 * Created on October 21, 2007, 2:15 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * a runtime exception that indicates a dimension mismatch
 * @author Sven Rahmann
 */
public class DimensionMismatchException extends java.lang.RuntimeException {
  
  /**
   * Creates a new instance of <code>DimensionMismatchException</code> without a detail message.
   */
  public DimensionMismatchException() {
  }
  
  
  /**
   * Constructs an instance of <code>DimensionMismatchException</code> with the specified detail message.
   * @param msg the detail message.
   */
  public DimensionMismatchException(String msg) {
    super(msg);
  }
}
