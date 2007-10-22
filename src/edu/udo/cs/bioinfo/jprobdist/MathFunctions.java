/*
 * MathFunctions.java
 *
 * Created on October 21, 2007, 3:32 AM
 *
 */

package edu.udo.cs.bioinfo.jprobdist;
import static java.lang.Math.*;

/**
 * This class contains only *static* Mathematical functions that
 * are useful for computing certain quantities related to probability distributions.
 * @author Sven Rahmann
 */
public final class MathFunctions {

  /** magic coefficients for computing the logarithmic gamma function */
  private final static double[] lngammacoeff =
  { 76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179E-02,
    -0.5395239384953E-05
  };
  
  /** Returns the logarithm of the Gamma function of the argument
   *  Throws an ArithmeticException if the argument is not positive.
   *@param xx  the argument of the logarithmic Gamma functon
   *@return the logarithm of the Gamma function of the argument
   */
  public static final
      double lngamma(final double xx) {
    double x,y,tmp,ser;
    if (xx<=0) throw new ArithmeticException("lngamma: Argument "+xx+" must be positive");
    y=x=xx;
    tmp=x+5.5;
    tmp-=(x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (int j=0; j<=5; j++) ser += lngammacoeff[j]/++y;
    return(-tmp+log(2.5066282746310005*ser/x));
  }
  
  /** returns the logarithm of x! (x factorial) */
  public static final
      double lnfactorial(final double x) {
    return(lngamma(x+1));
  }
  
  /** returns the logarithm of the binomial coefficient (n choose k) */
  public static final 
      double lnbincoeff(final double n, final double k) {
    return (lngamma(n+1) - lngamma(k+1) - lngamma(n-k+1));
  }
  
  /** returns (n choose k), the biomial coefficient of double arguments */
  public static final
      double bincoeff(final double n, final double k) {
    return exp(lnbincoeff(n,k));
  }
  
  /** returns (n choose k), the binomial coefficient of integer arguments */
  public static final
      long bincoeff(final int n, final int k) {
    return (long)xround(bincoeff(n,k));
  }
  
  /** rounds argument x to granularity eps */
  public static final
      double xround(final double x, final double eps) {
    final double y=eps*floor(x/eps+0.5);
    return (eps*floor(y/eps+0.5)); // doppelt haelt besser
  }
  
  /** rounds argument x to the nearest integer */
  public static final
      double xround(final double x) {
    return floor(x+0.5);
  }
  
  /** add two numbers in logarithmic space, i.e.,
   *  compute the logarithm of z = x+y, 
   *  given lx:=ln(x) and ly:=ln(y).
   *  Assuming lx>=ly, we use the formula 
   *  z = x + y  =  x*(1+y/x), i.e., 
   *  ln(z) = ln(x) + ln(1+y/x)
   *@param lx  the logarithm of x
   *@param ly  the logarithm of y
   *@return log(exp(lx)+epx(ly)) in a numerically stable way
   */
  public static final 
      double logsum(final double lx, final double ly) {
    if(lx>=ly) return (lx + log1p(exp(ly-lx)));
    else       return (ly + log1p(exp(lx-ly)));
  }
  
  
  /** computes n-th Fibonacci number F(n), where F(0)=F(1)=1,
   * and F(n) = F(n-1)+F(n-2) for n>=2.
   *@param n
   *@return F(n)
   */
  public static final long fib(int n) {
    if (n<0) throw new IllegalArgumentException("fib: parameter n must be nonnegative");
    if (n==0 || n==1) return 1L;
    long f2=1, f1=1, f=1;
    for (int i=2; i<=n; i++) {
      f = f1+f2;
      f2 = f1;
      f1 = f;
    }
    return f;
  }
  
  /** finds largest n such that the n-th Fibonacci number F(n)
   * is smaller or equal to a given number fmax
   * @param fmax  the upper bound on F(n)
   * @return largest n such that F(n)<=fmax
   */
  public static final int fibFind(long fmax) {
    if (fmax<=0) return -1;
    if (fmax==1) return 1;
    long f2=1, f1=1, f=2;
    int n;
    for (n=2; f<=fmax; n++) {
      f = f1+f2;
      f2 = f1;
      f1 = f;
    }
    return n-2;
  }
  
  /** computes n-th Fibonacci string FS(n), where FS(0)="a",
   *  FS(1)="b", and FS(n) = FS(n-1)+FS(n-2) for n>=2.
   *@param n
   *@return FS(n)
   */
  public static final String fibString(int n) {
    if (n<0) throw new IllegalArgumentException("fibString: parameter n must be nonnegative");
    if (n==0) return "b";
    if (n==1) return "a";
    String s2="b", s1="a", s=null;
    for (int i=2; i<=n; i++) {
      s = s1+s2;
      s2 = s1;
      s1 = s;
    }
    return s;
  }
  
// end class
}
