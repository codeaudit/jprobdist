/*
 * UnimodalDistribution.java
 *
 * Created on October 21, 2007, 6:20 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;

/**
 * a tagging interface that indicates that the implementing 
 * univariate (1-dimensional) distribution is unimodal, i.e.,
 * its pmf or density has only one local maximum (which may be 
 * an interval of consecutive values).
 *
 * @author Sven Rahmann
 */
public interface UnimodalDistribution extends UVDistribution {
  
  /** returns the location of the mode, 
   * or of the smallest mode in the case that the mode is an interval.
   */
  public double mode();
  
  /** returns the whole interval of modes */
  public Interval modeInterval();
  
}
