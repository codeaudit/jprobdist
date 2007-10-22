/*
 * BinomialDistributionTest.java
 * JUnit based test
 *
 * Created on October 21, 2007, 8:45 AM
 */

package edu.udo.cs.bioinfo.jprobdist;

import junit.framework.*;
import java.util.Iterator;
import static java.lang.Math.*;
import static edu.udo.cs.bioinfo.jprobdist.MathFunctions.*;
import java.util.NoSuchElementException;

/**
 *
 * @author rahmanns
 */
public class BinomialDistributionTest extends TestCase {
  
  public BinomialDistributionTest(String testName) {
    super(testName);
  }

  protected void setUp() throws Exception {
  }

  protected void tearDown() throws Exception {
  }

  public void testAll() {
    final int N = 10;
    BinomialDistribution b = new BinomialDistribution(N,0.4);
    double[] p = new double[N+1];
    for(int i=0; i<N; i++) {
      p[i] = b.cdf((double)i);
      System.out.printf("cdf(%d)=%f, qf(%f)=%f %n",i,p[i], p[i], b.qf(p[i]));
    }
  }


}
