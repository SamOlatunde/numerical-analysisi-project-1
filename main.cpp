// Huy Ngo, David Reyes
// Numerical Analysis
// Fall 2023
// The program finds the root of the equation that was set in the code
// using 4 different methods, Newton's, Bisection, Fixed Point, and Secant
// it also tell you the solution after X iterations or failure.

#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
using namespace std;

const double PI = 3.14;

struct RootFinding {
  function<double(double)> f;  // The function for which to find the root
  function<double(double)> df; // The derivative of the function
  function<double(double)> g;

  RootFinding(function<double(double)> func, function<double(double)> dfunc,
              function<double(double)> gfunc)
      : f(func), df(dfunc), g(gfunc) {}

  // BISECTION METHOD
  // INPUT endpoints a, b; tolerance TOL; maximum number of iterations N0.
  // OUTPUT approximate solution p or message of failure.
  // Step 1  Set i = 1;
  //         FA = f (a).
  // Step 2  While i ≤ N0 do Steps 3–6.
  // Step 3  Set p = a + (b − a)/2; (Compute pi.)
  //         FP = f ( p).
  // Step 4  If FP = 0 or (b − a)/2 < TOL then
  //         OUTPUT (p); (Procedure completed successfully.)
  //         STOP.
  // Step 5  Set i = i + 1.
  // Step 6  If FA · FP > 0 then set a = p; (Compute ai, bi.)
  //         FA = FP
  //         else set b = p. (FA is unchanged.)
  // Step 7  OUTPUT (‘Method failed after N0 iterations, N0 =’, N0);
  //         (The procedure was unsuccessful.)
  //         STOP
  double bisectionMethod(double a, double b, double tol, int maxIter) {
    // Step 1
    int i = 1;
    double FA = f(a);
    double p, FP;

    // Step 2
    while (i <= maxIter) {

      // Step 3
      p = a + (b - a) / 2;
      FP = f(p);

      // Step 4
      if (FP == 0.0 || (b - a) / 2 < tol) {
        cout << "Bisection Method approximate solution after " << i
             << " iterations: " << p << endl;
        return p;
      }

      // Step 5 and Step 6
      if (FA * FP > 0) {
        a = p;
        FA = FP;
      } else {
        b = p;
      }

      i++;
    }

    // Step 7
    cout << "Bisection method failed after " << maxIter << " iterations."
         << endl;
    return p;
  }

  // NEWTON'S METHOD
  // INPUT  initial approximation p0; tolerance TOL; maximum number of
  // iterations N0.
  // OUTPUT  approximate solution p or message of failure.
  // Step 1  Set i = 1.
  // Step 2  While i ≤ N0 do Steps 3–6.
  // Step 3  Set p = p0 − f (p0)/f ( p0). (Compute pi.)
  // Step 4  If | p − p0| < TOL then
  //         OUTPUT (p); (The procedure was successful.)
  //         STOP.

  // Step 5  Set i = i + 1.
  // Step 6  Set p0 = p. (Update p0.)
  // Step 7  OUTPUT (‘The method failed after N0 iterations, N0 =’, N0);
  //         (The procedure was unsuccessful.)
  //         STOP.
  double newtonMethod(double p0, double tol, int maxIter) {
    // Step 1
    int i = 1;
    double p;
    // Step 2
    while (i <= maxIter) {
      // Step 3
      p = p0 - f(p0) / df(p0);
      // Step 4
      if (fabs(p - p0) < tol) {
        cout << "Newton's method approximate solution after " << i
             << " iterations: " << p << endl;
        return p;
      }

      i++;    // Step 5
      p0 = p; // Step 6
    }

    cout << "Newton's method failed after " << maxIter << " iterations."
         << endl; // Step 7
    return p0;
  }

  // SECANT METHOD
  // INPUT   initial approximations p0, p1; tolerance TOL; maximum number of
  // iterations N0.
  // OUTPUT  approximate solution p or message of failure.
  // Step 1  Set i = 2; q0 = f ( p0); q1 = f ( p1).
  // Step 2  While i ≤ N0 do Steps 3–6.
  // Step 3  Set p = p1 − q1( p1 − p0)/(q1 − q0). (Compute pi.)
  // Step 4  If | p − p1| < TOL then
  //         OUTPUT (p); (The procedure was successful.)
  //         STOP.
  // Step 5  Set i = i + 1.
  // Step 6  Set p0 = p1; (Update p0, q0, p1, q1.)
  //         q0 = q1;
  //         p1 = p;
  //         q1 = f ( p).
  // Step 7  OUTPUT (‘The method failed after N0 iterations, N0 =’, N0);
  //         (The procedure was unsuccessful.)
  //         STOP.
  double secantMethod(double p0, double p1, double tol, int maxIter) {
    // Step 1
    int i = 2;
    double q0 = f(p0);
    double q1 = f(p1);

    // Step 2
    while (i <= maxIter) {
      // Step 3
      double p = p1 - q1 * (p1 - p0) / (q1 - q0);

      // Step 4
      if (fabs(p - p1) < tol) {
        cout << "Secant method approximate solution after " << i - 1
             << " iterations: " << p << endl;
        return p;
      }

      // Step 5
      i++;

      // Step 6
      p0 = p1;
      q0 = q1;
      p1 = p;
      q1 = f(p);
    }

    // Step 7
    cout << "Secant method failed after " << maxIter << " iterations." << endl;
    return p1;
  }

  // FIXED POINT METHOD
  // INPUT   initial approximation p0; tolerance TOL; maximum number of
  // iterations N0.
  // OUTPUT  approximate solution p or message of failure.
  // Step 1  Set i = 1.
  // Step 2  While i ≤ N0 do Steps 3–6. Step 3  Set p = g( p0).
  //         (Compute pi.)
  // Step 4  If | p − p0| < TOL then
  //           OUTPUT ( p); (The procedure was successful.)
  //           STOP.
  // Step 5  Set i = i + 1.
  // Step 6  Set p0 = p. (Update p0.)
  double fixedPointMethod(double p0, double tol, int maxIter) {
    // Step 1
    int i = 1;

    // Step 2
    while (i <= maxIter) {
      // Step 3
      double p = g(p0);

      // Step 4
      if (fabs(p - p0) < tol) {
        cout << "Fixed point found after " << i << " iterations: " << p << endl;
        return p;
      }

      // Step 5 and Step 6
      i++;
      p0 = p;
    }

    cout << "Max iterations reached! Fixed point may not have been found."
         << endl;
    return p0;
  }
};

void testRootMethods(RootFinding &rf, const string &title, double init1,
                     double init2, double tol, int maxIter) {
  cout << "Testing function: " << title << endl;

  auto start = std::chrono::high_resolution_clock::now();
  rf.newtonMethod(init1, tol, maxIter);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  // Format the time taken to display it in .0000000 format
  std::cout << std::fixed << std::setprecision(7);

  cout << "Time taken by Newton's Method: " << elapsed.count() << " seconds"
       << endl
       << endl;

  start = std::chrono::high_resolution_clock::now();
  rf.bisectionMethod(init1, init2, tol, maxIter);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  cout << "Time taken by Bisection Method: " << elapsed.count() << " seconds"
       << endl
       << endl;

  start = std::chrono::high_resolution_clock::now();
  rf.secantMethod(init1, init2, tol, maxIter);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  cout << "Time taken by Secant Method: " << elapsed.count() << " seconds"
       << endl
       << endl;

  start = std::chrono::high_resolution_clock::now();
  rf.fixedPointMethod(init1, tol, maxIter);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  cout << "Time taken by Fixed Point Method: " << elapsed.count() << " seconds"
       << endl
       << endl;

  cout << "-------------------------------------------------" << endl;
}

int main() {
  double tol = 1e-6;
  int maxIter = 100;

  // Example 1: f(x) = x^2 - 6
  // set g(x) for fixed point method = Sqrt(6)
  RootFinding example1([](double x) { return x * x - 6; },
                       [](double x) { return 2 * x; },
                       [](double x) { return sqrt(6); });
  testRootMethods(example1, "x^2 - 6", 1.0, 3.0, tol, maxIter);

  // Example 2: f(x) = e^x - 3
  // set g(x) for fixed point method = 1 - PI
  RootFinding example2([](double x) { return exp(x) - 3; },
                       [](double x) { return exp(x); },
                       [](double x) { return log(3); });
  testRootMethods(example2, "e^x - 3", 1.0, 2.0, tol, maxIter);

  // Example 3: f(x) = e^x - 3x^2
  // set g(x) for fixed point method = ln(3x^2)
  RootFinding example3([](double x) { return exp(x) - 3 * x * x; },
                       [](double x) { return exp(x) - 6 * x; },
                       [](double x) { return log(3 * pow(x, 2)); });
  testRootMethods(example3, "e^x - 3x^2", 0.5, 3.0, tol, maxIter);

  return 0;
}
