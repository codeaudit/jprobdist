/*
 * MathFunctionsTest.java
 * JUnit based test
 *
 * Created on October 22, 2007, 7:28 PM
 */

package edu.udo.cs.bioinfo.jprobdist;

import junit.framework.*;
import static java.lang.Math.*;

/**
 *
 * @author rahmanns
 */
public class MathFunctionsTest extends TestCase {
  
  public MathFunctionsTest(String testName) {
    super(testName);
  }

  protected void setUp() throws Exception {
  }

  protected void tearDown() throws Exception {
  }

  /**
   * Test of lngamma method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testLngamma() {
    System.out.println("lngamma");
    
    double xx = 0.5;
    
    double expResult = Math.log(Math.PI)/2.0;
    double result = MathFunctions.lngamma(xx);
    assertTrue(String.format("Expected %.16f, Result=%.16f",expResult,result), 
        MathFunctions.DoubleEqual(expResult, result));   
  }

  /**
   * Test of lnfactorial method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testLnfactorial() {
    System.out.println("lnfactorial");
    
    double x = 4.0;
    
    double expResult = Math.log(24);
    double result = MathFunctions.lnfactorial(x);
    assertTrue(String.format("Expected %.16f, Result=%.16f",expResult,result), 
        MathFunctions.DoubleEqual(expResult, result));   
  }

  /**
   * Test of lnbincoeff method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testLnbincoeff() {
    System.out.println("lnbincoeff");
    
    double n = 10;
    double k = 2;
    
    double expResult = Math.log(45);
    double result = MathFunctions.lnbincoeff(n, k);
    assertTrue(String.format("Expected %.16f, Result=%.16f",expResult,result), 
        MathFunctions.DoubleEqual(expResult, result));   
  }

  /**
   * Test of bincoeff method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testBincoeff() {
    System.out.println("bincoeff");
    
    int n = 10;
    int k = 2;
    
    long expResult = 45;
    long result = MathFunctions.bincoeffI(n, k);
    assertEquals(expResult, result);
  }

  /**
   * Test of xround method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testXround() {
    System.out.println("xround");
    
    double x = -0.6;
    double eps = 1.0;
    
    double expResult = -1.0;
    double result = MathFunctions.xround(x, eps);
    assertEquals(expResult, result);
  }

  /**
   * Test of logsum method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testLogsum() {
    System.out.println("logsum");
    
    double lx = Double.NEGATIVE_INFINITY; // 0.0
    double ly = -0.5;                     // exp(-0.5)
    double expResult = ly;
    double result = MathFunctions.logsum(lx, ly);
    assertEquals(expResult, result);
    
    // Try 0+0=0
    lx = Double.NEGATIVE_INFINITY; 
    ly = lx;                    
    expResult = ly;
    result = MathFunctions.logsum(lx, ly);
    assertEquals(expResult, result);
    
  }


  /**
   * Test of findRootInterval method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testFindRootInterval() {
    System.out.println("findRootInterval");
    
    MathFunctions.RealFunction f = null;
    double a = 0.0;
    double b = 0.0;
    
    Interval expResult = null;
    Interval result = MathFunctions.findRootInterval(f, a, b);
    assertEquals(expResult, result);
    
    // TODO review the generated test code and remove the default call to fail.
    fail("The test case is a prototype.");
  }

  /**
   * Test of findRootBisection method, of class edu.udo.cs.bioinfo.jprobdist.MathFunctions.
   */
  public void testFindRootBisection() {
    System.out.println("findRootBisection");
    
    MathFunctions.RealFunction f = null;
    double a = 0.0;
    double b = 0.0;
    double xacc = 0.0;
    
    double expResult = 0.0;
    double result = MathFunctions.findRootBisection(f, a, b, xacc);
    assertEquals(expResult, result);
    
    // TODO review the generated test code and remove the default call to fail.
    fail("The test case is a prototype.");
  }

  
}
