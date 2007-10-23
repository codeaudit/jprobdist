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


/**
 *
 * @author Sven Rahmann
 */
public class BinomialDistributionTest extends TestCase {
  
  public BinomialDistributionTest(String testName) {
    super(testName);
  }

  protected void setUp() throws Exception {
  }

  protected void tearDown() throws Exception {
  }

  public void testSmoke() {
    BinomialDistribution b = null;
    final int N = 20;
    
    // test cdf() in boundary cases
    for (double p=0; p<=1; p+=1) {
      b = new BinomialDistribution(N,p);
      for(int i=0; i<=N; i++)  assertEquals(i==N*p? 1.0: 0.0, b.P(i));
      for(int i=0; i<=N; i++)  assertEquals(i>=N*p? 1.0: 0.0, b.cdf(i));
    }
    
    
    // test qf()
    b = new BinomialDistribution(N,0.4);
    for(int i=0; i<=N; i++) {
      double pp = b.cdf((double)i);
      double qq = b.qf(pp);
      assertEquals("",(double)i,qq,0.0);
    }
    for(double p=0.0; p<=1; p+=1.0/128) {
      double qq = b.qf(p);
      double pp = b.cdf(qq);
      double ppp = b.cdf(qq-1);
      //System.out.printf("p=%f, qf(p)=%f, cdf(qf(p))=%f, cdf(qf(p)-1)=%f%n", p,qq,pp,ppp);
      assertTrue("",p<=pp);
      assertTrue("",p>ppp || (p==0.0 && p>=ppp));
    }
    
    // test median vs qf
    for(double p=0.0; p<=1; p+=1.0/128) {
      b = new BinomialDistribution(N,p);
      double m = b.median();
      double q = b.qf(0.5);
      //System.out.printf("median = %f and %f%n", m, q);
      assertEquals("", m,q, 0.0);   
    }
    
  }
  
  

}
