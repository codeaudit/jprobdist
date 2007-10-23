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

  
  //============ GAMMA, BETA, FACTORIAL, BINOMIAL COEFF, etc. ================
  
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
    if (xx<=0) throw new ArithmeticException("lngamma: Argument "+xx+" > 0 required");
    y=x=xx;
    tmp = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015; // nine zeros
    for (int j=0; j<=5; j++) ser += lngammacoeff[j]/++y;
    return(-tmp+log(2.5066282746310005*ser/x));
  }
  
  private static final double[] factorialTable   = new double[33];
  private static final double[] lnfactorialTable = new double[101];
  private static int factorialTableFill = -1;
  
  /** returns n!, i.e., n factorial, as a double number */
  public static final synchronized
      double factorial(final int n) {
    if(n<0) throw new IllegalArgumentException("factorial: Argument "+n+" >= 0 required");
    if(n>32) return(xround(exp(lngamma(n+1.0))));
    // n is an integer in {0, 1, ..., 32}
    if (factorialTableFill==-1) { factorialTable[0]=factorialTable[1]=1.0; factorialTableFill=1; }
    while(factorialTableFill < n) {
      final int j = factorialTableFill++;
      factorialTable[factorialTableFill] = factorialTable[j]*factorialTableFill;
    }
    return factorialTable[n];
  }
  
  
  /** returns ln(x!), i.e. the logarithm of x factorial, for x &gt; -1 */
  public static final
      double lnfactorial(final double x) {
    if(x<=-1.0) throw new IllegalArgumentException("lnfactorial: Argument "+x+" > -1 required");
    if(x!=round(x) || x>100) return(lngamma(x+1));
    // x is an integer in {0, 1, ..., 100}
    if(x<=1.0) return 0.0;
    final int xi = (int)x;
    return (lnfactorialTable[xi]!=0.0)? lnfactorialTable[xi] 
        : (lnfactorialTable[xi] = lngamma(x+1));
  }
  
  /** returns the logarithm of the binomial coefficient (n choose k) */
  public static final 
      double lnbincoeff(final double n, final double k) {
    return (lnfactorial(n) - lnfactorial(k) - lnfactorial(n-k));
  }
  
  /** returns (n choose k), the biomial coefficient,
   *  of the given double arguments
   */
  public static final
      double bincoeff(final double n, final double k) {
    return exp(lnbincoeff(n,k));
  }
 
  /** returns (n choose k), the biomial coefficient,
   *  of the given double arguments,
   *  rounded to the next integer, as a double
   */
  public static final
      double bincoeffR(final double n, final double k) {
    return xround(exp(lnbincoeff(n,k)));
  }
  
  /** returns (n choose k), the binomial coefficient of integer arguments,
   *  as a long integer. In contrast to bincoeff(n,k), this function provides
   *  integer rounding at the end, but it also requires integral arguments.
   */
  public static final
      long bincoeffI(final int n, final int k) {
    //TODO: Deal with negative n, k -- see Concrete Mathematics
    return round(exp(lnbincoeff(n,k))); // autoconversion to double for n,k
  }
  
  /** returns the incomplete Gamma function P(a,x).
   *  This is the cdf of a Gamma distribution with shape a and scale 1.
   */
  public static final 
      double gammaP(final double a, final double x) {
    if (x<0.0 || a<=0.0) throw new IllegalArgumentException("a>0 and x>=0 required");
    if (x < a+1.0) return gammaSeries(a,x)[0];
    else           return 1.0 - gammaCF(a,x)[0];
  }

  /** returns the incomplete Gamma function Q(a,x).
   *  This is the ucdf of a Gamma distribution with shape a and scale 1.
   */
  public static final 
      double gammaQ(final double a, final double x) {
    if (x<0.0 || a<=0.0) throw new IllegalArgumentException("a>0 and x>=0 required");
    if (x < a+1.0) return 1.0 - gammaSeries(a,x)[0];
    else           return gammaCF(a,x)[0];
  }
  
  /** computes the incomplete gamma function gammaP(a,x) via series representation.
   *  via continued fractions representation.
   * Use this function when x <= a+1.
   *@param a the 'shape parameter'
   *@param x
   *@return a double array gammaSeries[] with gammaSeries[0] = gammaP(a,x) and
   * gammaSeries[1] = lngamma(a)   
   */
  public static final
      double[] gammaSeries(final double a, final double x) {
    final double ITMAX = 200; // max number of iterations
    final double gln = lngamma(a);
    if (x<=0.0) {
      if (x==0.0) return new double[] {0.0, gln};
      throw new IllegalArgumentException("x >=0 required");
    }
    // now x>0, a>0; otherwise we would have been kicked out
    double ap = a;
    double del = 1.0/a;
    double sum = del;
    for(int n=1; n<=ITMAX; n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if(abs(del)<abs(sum)*DBL_TOL) {
        return new double[] {sum*exp(-x+a*log(x)-gln), gln};
      }
    }
    throw new RuntimeException("a too large, ITMAX too small");
  }

  /** computes the upper incomplete gamma function gammaQ(a,x) 
   *  via continued fractions representation.
   * Use this function when x >= a+1.
   *@param a the 'shape parameter'
   *@param x
   *@return a double array gammaCF with gammaCF[0]==gammaQ(a,x) and
   * gammaCF[1] = lngamma(a)
   */
  public static final
      double[] gammaCF(final double a, final double x) {
    final double ITMAX = 200; // max number of iterations
    final double gln = lngamma(a);
    double b = x+1.0-a;
    double c = 1.0/DBL_MIN_NORMAL;
    double d = 1.0/b;
    double h = d;
    int i;
    for(i=0; i<=ITMAX; i++) {
      final double an = -i*(i-a);
      b += 2.0;
      d = an*d+b;  if (abs(d)<DBL_MIN_NORMAL) d = DBL_MIN_NORMAL;
      c = b+an/c;  if (abs(c)<DBL_MIN_NORMAL) c = DBL_MIN_NORMAL;
      d = 1.0/d;
      final double del = d*c;
      h *= del;
      if (abs(del-1.0)<2*DBL_TOL) break;
    }
    if (i>ITMAX) throw new RuntimeException("a too large; ITMAX too small");
    return new double[] { exp(-x+a*log(x)-gln)*h, gln };
  }

  
  /** returns the logarithm of the beta function at its arguments */
  public static final 
      double lnbeta(final double z, final double w) {
    return (lngamma(z)+lngamma(w)-lngamma(z+w));
  }
  
  /** returns the beta function at its arguments */
  public static final 
      double beta(final double z, final double w) {
    return exp(lngamma(z)+lngamma(w)-lngamma(z+w));
  }
  
  
  
  
  //================= Round to Integer, return as double ===================
  
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
  
  
  //================== Logspace arithmetic ==================================
  
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
    if(lx>=ly) {
      if (ly==Double.NEGATIVE_INFINITY) return lx;
      return (lx + log1p(exp(ly-lx)));
    }
    else {
      if (lx==Double.NEGATIVE_INFINITY) return ly;
      return (ly + log1p(exp(lx-ly)));
    }
  }
  
  /** add numbers in logarithmic space, i.e.,
   *  compute the logarithm of sum(v1,v2,v3,...),
   *  when ln(v1), ln(v2), ln(v3), ... are given
   */
  public static final double logsum(final double... lv) {
    if (lv.length<=2) { 
      if (lv.length==2) return logsum(lv[0],lv[1]);
      if (lv.length==1) return lv[0];
      if (lv.length==0) return Double.NEGATIVE_INFINITY;
    }
    // at least 3 values, find largest one, call it m, and its position mj
    double m  = Double.NEGATIVE_INFINITY;
    int    mj = 0;
    for(int j=0; j<lv.length; j++) 
      if(lv[j]>m) { m=lv[j]; mj=j; }
    if (m==Double.NEGATIVE_INFINITY) return m; // all lv[j] are -inf
    // sum up  exp(lv[j]-m) for j != mj
    double sum=0;
    for(int j=0;    j<mj;        j++) sum+=exp(lv[j]-m);
    for(int j=mj+1; j<lv.length; j++) sum+=exp(lv[j]-m);
    return ( m + log1p(sum));
  }
      
  
  //================== FIBONACCI NUMBERS AND STRINGS =========================
  
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
    //TODO: Use a string buffer instead to do this more efficiently
    for (int i=2; i<=n; i++) {
      s = s1+s2;
      s2 = s1;
      s1 = s;
    }
    return s;
  }
  
  //================= REAL FUNCTIONS ===========================================
  
  /** this interface defines a real function x->f(x) */
  public interface RealFunction {
    /** a real function that maps x to f(x) */
    public double valueAt(double x);
  }
  
  /** the identity function */
  public static final RealFunction ID = new RealFunction() {
    public double valueAt(double x) {return x;}
  };

  
  //================= ROOT FINDING ===========================================
  /** find an interval that contains a root of f, starting with (a,b) and 
   *expanding to the left and right as necessary */
  public static final Interval findRootInterval(final RealFunction f, double a, double b) {
    if (a>=b || Double.isNaN(a) || Double.isNaN(b)) throw new IllegalArgumentException("requires a<b");
    if (a<=Double.NEGATIVE_INFINITY) a = -(0x1.fffffffffffffP+511);
    if (b>=Double.POSITIVE_INFINITY) b =   0x1.fffffffffffffP+511;
    double fa = f.valueAt(a);
    double fb = f.valueAt(b);
    for(int j=0; j<50; j++) {
      if (fa*fb<0.0) return new Interval(a,b,Interval.Type.Open);
      if (abs(fa)<abs(fb)) fa = f.valueAt(a += 1.6*(a-b));
      else                 fb = f.valueAt(b += 1.6*(b-a));
    }
    return Interval.TheNaNInterval;
  }
  
  
  
  /** find a root of f in the interval [a,b] using bisection.
   * If there is a single interval of roots, find the smallest root.
   *@param f  the function
   *@param a  left boundary of interval
   *@param b  right boundary of interval
   *@param xacc desired x-accuracy: the function returns as soon as the length
   * of the interval drops below xacc. Specify xacc=0.0 for maximum precision!
   */
  public static final double findRootBisection(final RealFunction f, double a, double b, final double xacc) {
    double fa,fb;
    
    fa = f.valueAt(a); if(fa==0.0) return a;
    fb = f.valueAt(b);
    if (fa*fb>0.0) throw new IllegalArgumentException(
        String.format("f(a), f(b) must have opposite signs: [a,b]=[%f,%f], f(a)=%f, f(b)=%f", a,b,fa,fb));
    
    double x, fx, dx = b-a; 
    do {
      fx = f.valueAt(x=(a+(dx*=0.5)));
      if(dx<xacc) return x;
      if (fx<0) {a=x;} else {b=x;}
    } while (!DoubleEqual(a,b));
    return (a+b)/2.0;
  }
  
  
  public final static double DBL_EPS = 0x1.0p-52;
  public final static double TOL_FACTOR = 2;
  public final static double DBL_TOL = TOL_FACTOR*DBL_EPS;
  public final static double DBL_MIN_NORMAL = 0x1.0p-1022;

  
  /** returns min {|x| : x in [a,b]} * DBL_TOL  */
  public static final double DoubleAccuracyIn(double a, double b) {
    if(b<a) {final double c=a; a=b; b=c;}  // ensures that now a<=b
    final double m = (b<0)? -b : ( (a>=0)? a : ((-a<b)? -a:b) ); // m=min(|a|,|b|)
    return (m<=DBL_MIN_NORMAL)? TOL_FACTOR*DBL_MIN_NORMAL : DBL_TOL*m;
  }
  
  /** checks for two double numbers a<b if numerically a==b, 
   * i.e., if (b-a) <= DBL_TOL*max(|a|,|b|)
   * This REQUIRES that a<=b, otherwise an incorrect value may result.
   */
  public static final boolean DoubleEqual(final double a, final double b) {
    //if(b<a) {final double c=a; a=b; b=c;} // now a<=b // we skip this!
    final double m = (b<0)? -a : ( (a>=0)? b : ((-a>b)? -a:b) ); // m=max(|a|,|b|)
    return (m<= TOL_FACTOR*DBL_MIN_NORMAL)? (b-a <= TOL_FACTOR*Double.MIN_VALUE) : (b-a <= DBL_TOL*m);
  }
  
// end class
}
