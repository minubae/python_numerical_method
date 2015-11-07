# Numerical Method (Analysis) Module - pynumerical.py
# Date: Oct/21/2015, Wed - 
# Version: ver 0.1
# Author: Minwoo Bae
# Contact: minubae.nyc@gmail.com

import math
from sympy import *

## A Function for the purpose of testing
fx = lambda x: 2**x - 3
fx2 = lambda x:  x**3 + 4*x**2 - 10
fx3 = lambda t,y: y/t - (y/t)**2  # --> Differential Equation; for testing of this Euler's method.
fx4 = lambda x: cos(x)-x

### 1. Mathematical Preliminaries and Error Analysis

## Decimal Machine Numbers: The actual error; the absolute error; the relative error
# Suppose that a_p is an approximation to p.
# The actual error is p-a_p,
# the abosulte error is |p-a_p|,
# and the relative error is |p-a_p|/|p|,provided that p != 0.
def absolute_value(p):
        return math.fabs(p)

def actual_error(p, a_p):
        return p - a_p

def absolute_error(p, a_p):
        return math.fabs(p-a_p)

def relative_error(p, a_p):
        result = 0
        if absolute_value(p) != 0:
                result = absolute_error(p,a_p)/absolute_value(p)
                return result
        else:
                print('error: absolute_value must not be equal to ',absolute_value(p))

### 2. Solution of Equations in One Variable
                
## The Bisection Method (Binary Search)
# The Bisection, or Binary-search, Method is based on the Intermediate Value Theorem.
# Suppose f is a continuous function defined on the interval [a,b], with f(a) and f(b)
# of opposite sign. THe Intermediate Value Theorem implies that a number p exists in (a,b)
# with f(p) = 0. Although the procedure will work when there is more than one root in the
# interval (a,b), we assume for simplicity that the root in this interval is unique.
# The method calls for a repeated halving (or bisecting) of subintervals of [a,b] and,
# at each step, locating the half containing p.
# To begin, set a1 = a and b1 = b and let p1 be the midpoint of [a,b]; that is,
# p1 = a1 + (b1-a1)/2 = (a1+b1)/2.
# If f(p1) = 0, then p = p1, and we are done
# if f(p1) != 0, then f(p1) has the same sign as either f(a1) or f(b1).
# --- if f(p1) and f(a1) have the same sign, p is an element of (p1,b1). Set a2=p1 and b2=b1.
# --- if f(p1) and f(a1) have opposite sign, p is an element of (a1,p1). Set a2=p1 and b2=p1.
# INPUT: Function fx; endpoints a, b; tolerance TOL; maximum number of iterations N0
# OUTPUT: Approximate solution p or message of failure.
def bisection_method(fx, a, b, tol, N):
        i = 1
        FA = fx(a)

        try:
            while i <= N:

                p = a + (b-a)/2
                FP = fx(p)

                print(i,': ', 'p:', p, '; f(p):', FP)
                
                if FP == 0 or (b-a)/2 < tol:
                    return p; break
                i+=1
                
                if FA*FP > 0:
                        a = p
                else:
                        b = p
        except:
            return 'The procedure was unsuccessful.'
                
# def bisection_method(fx, a, b, num):
#         i = 0
#         pivot = 0
#         isNegative = lambda fx, a, b: True if fx(a) * fx(b) < 0 else False
#         findPivot = lambda a, b: (a+b)/2
        
#         if isNegative(fx, a, b):
#                 for i in range(num):
#                         pivot = findPivot(a, b)
#                         if fx(pivot) > 0 and fx(b) < 0:
#                                 a = pivot
#                         else:
#                                 b = pivot
#                 return pivot
#         else:
#                 print('Please reset up the right interval.')

## Fixed-Point Iteration Theorem
# Let g is continuous on [a,b] be such that g(x) exists on [a,b], for all x in [a,b].
# Suppose, in addition, that g' exists on (a,b) and that a constant 0 < k < 1 exists with
# |g'(x)| <= k, for all x in (a,b).
# Then, for any number p0 in [a,b], the sequence defiend by p_n = g(p_n-1), n >= 1, 
# converges to the unique fixed point p in [a,b].
# INPUT: Initial approximation p0; tolerance TOL; maximum number of iteration N.
# OUTPUT: Approximate solution p or message of failure.   
def fixed_point(p0, g, tol, N):
    i = 1
    try:
        while i <= N:
            p = g(p0)

            print('i:',i, 'p:',p)

            if math.fabs(p-p0) < tol:
                return p; break

            i += 1
            p0 = p

    except:
        return 'The procedure was unsuccessful'

## Newton's (or the Newton-Raphson)Method
# Newton's (or the Newton-Raphson)Method is one of the most powerful and well-known numerical methods for solving a root-finding problem.
# Newton's method is a functional iteration technique with p_n, for which
# P_n = P_n-1 - f(P_n-1)/f'(P_n-1), for n >= 1
# INPUT: Initial approximation p0; tolerance TOL; maximum number of iterations N.
# OUTPUT: Approximation solution p or message of failure.
x = symbols('x')
f = lambda x: math.cos(x)-x

def newton_method(p0, f, tol, N):
        i = 1
        try:
                while i <= N:
                        # TODO:
                        # diff(f(p0), x) needs to be fixed.
                        p = p0 - f(p0)/diff(f(p0),x)

                        if math.fabs(p-p0) < tol:
                                return p; break
                        i+=1
                        p0 = p
        except:
                return 'The procedure was unsuccessful'

## The Secant Method
# Newton's method is an extremly powerful technique, but it has a major weakness: the need to know
# the value of the derivative of f at each approximation. Frequently, f'(x) is far more difficult and needs
# more arithmetic operations to calculate than f(x).
# P_n = P_n-1 - (f(P_n-1)(P_n-1 - P_n-2)) / (f(P_n-1) - f(P_n-2))
# INPUT: Initial approximation p0, p1; tolerance TOL; maximum number of iterations N.
# OUTPUT: Approximate solution p or message of failure.
def secant_method(p0,p1,f,tol,N):
        i = 2
        q0 = f(p0)
        q1 = f(p1)
        
        try:
            while i <= N:
                p = p1 - q1*(p1-p0) / (q1-q0)
                print('i:',i, 'p:',p)

                if math.fabs(p-p1)< tol:
                    return p; break
                i+=1
                        
                p0 = p1
                q0 = q1
                p1 = p
                q1 = f(p)

        except:
                return 'The procedure was unsuccessful'

                        
## The Method of False Postion
# The Method of False Position (also called Regula Falsi) generates approximations in the same manner as the Secant
# method, but it includes a test to ensure that the root is always bracketed between successive iterations.
# P_n = P_n-1 - (f(P_n-1)(P_n-1 - P_n-2)) / (f(P_n-1) - f(P_n-2))
# INPUT: Initial approximation p0, p1; tolerance TOL; maximum number of iterations N.
# OUTPUT: Approximate solution p or message of failure.
def false_position(p0, p1, f, tol, N):

        i = 2
        q0 = f(p0)
        q1 = f(p1)

        try:
                while i <= N:
                        p = p1 - q1*(p1-p0)/(q1 - q0)
                        print('i:',i, 'p:',p)
                        
                        if math.fabs(p - p1) < tol:
                                return p; break
                        i+=1
                        q = f(p)

                        if q*q1 < 0:
                                p0 = p1
                                q0 = q1
                        p1 = p
                        q1 = q
        except:
                return 'The procedure was unsuccessful'

### 3. Interpolation and Polynomial Approximation
        
## Neville's Iterated Interpolation
# Theorem 3.5: Let f be defined at x0, x1,...,xk and let xj and xi be two distinct numbers in this set.
# Then P(x) = (x-xj)P0,1,...,j-1,j+1,...,k(x) - (x-xi)P0,1,...,j-1,j+1,...,k(x) / (xi-xj)
# is the kth Lagrange polynomial that interpolates f at the k+1 points x0, x1,..., xk.
# The procedure that uses the result of Theorem 3.5 to recursively generate interpolating
# polynomial approximations is called Nevill's method.
# Let Qi,j(x), for 0 <= j <= i, denote the interpolating polynomial of degree j on the (j+1) numbers
# xi-j, xi-j+1,..., xi-1,xi; that is, Qi,j = Pi-j, i-j+1,.., i-1, i.

# To evaluate the interpolating polynomial P on the n+1 distinct numbers x0,...xn
# at the number x for the function f
# INPUT: Numbers x, x0, x1,...,xn; values f(x0), f(x1),...,f(xn) as the first column Q0,0, Q1,0,...,Qn,0 of Q.
# OUTPUT: The table Q with P(x) = Qn,n.
def neville_method(x0, y, fx):

        n = len(x)        
        Q = [[0 for i in range(n)] for j in range(n)]

        for i in range(n):
                Q[i][0] = fx[i]
                
        for i in range(1,n):
                for j in range(1, i+1):
                        Q[i][j] = ((x0-x[i-j])*Q[i][j-1] - (x0 -x[i])*Q[i-1][j-1]) / (x[i] - x[i-j])

        print("Result of Neville's Method: \n")
        for i in range(n):
                for j in range(n):
                        print("%.7f" %Q[i][j])
                print("\n")

        return Q

## Newton's Divided-Difference Formula
# Iterated interpolation was used in the previous section to generate successively higher-degree
# polynomial approximations at a specific point. Divided-difference methods are used to successively
# generate the polynomial themselves.
# To obtain the divided-difference coefficients of the interpolatory polynomial P on the (n+1)
# distinct numbers x0, x1, ... , xn for the function f:
# INPUT: Numbers x0, x1, ... , xn; values f(x0), f(x1),...,f(xn) as F0,0, F1,0,...,Fn,0.
# OUTPUT: The numbers F0,0, F1,1,...,Fn,n where
# Pn(x) = F0,0 + Sum (Fi,i) from i=1 to n * Product(x-xj) from j=0 to i-1. (Fi,i is f[x0,x1,...,xi])
def divided_differences(x, fx):
        # x = [1.0, 1.3, 1.6, 1.9, 2.2]
        # fx = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
        n = len(x)
        F = [ [ 0 for i in range(n) ] for j in range(n) ]
        
        for i in range(n):
                F[i][0] = fx[i]
        
        for i in range(1, n):
                for j in range(1, i+1):
                        F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (x[i] - x[i-j])
                        
        return F

## Hermote Interpolation
# To obtain the coefficients of the Hermite interpolating polynomial H(x) on the (n+1) distinct numbers x0,..,xn
# for the function f :
# INPUT: Numbers x0, x1, ... , xn; values f(x0),...,f(xn) and f'(x0),...,f'(xn).
# OUTPUT: The numbers Q0,0, Q1,1,...,Q2n+1,2n+1 where
# H(x) = Q0,0 + Q1,1(x-x0) + Q2,2(x-x0)^2 + Q3,3(x-x0)^2(x-x1) + Q4,4(x-x0)^2(x-x1)^2 + ...
#          + Q2n+1,2n+1(x-x0)^2(x-x1)^2 ... (x-xn-1)^2(x-xn).
def hermite_interpolation(x, fx, fp):
        #x = [1.3, 1.6, 1.9]
        #fx = [0.6200860, 0.4554022, 0.2818186]
        #fp = [-0.5220232, -0.5698959, -0.5811571]
        
        n = len(x)
        z = [0 for i in range(2*n)]
        Q = [ [ 0 for i in range(2*n) ] for j in range(2*n) ]
        #z = [0]*(2*n) #Q = [[None]*(2*n)]*(2*n)

        for i in range(n):
                z[2*i] = x[i]
                z[2*i+1] = x[i]
                Q[2*i][0] = fx[i]
                Q[2*i+1][0] = fx[i]
                Q[2*i+1][1] = fp[i]
                if i != 0:
                        Q[2*i][1] = (Q[2*i][0] - Q[2*i-1][0]) / (z[2*i] - z[2*i-1])

        for i in range(2, 2*n):
                for j in range(2, i+1):
                        Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (z[i] - z[i-j])
                        
        return Q

### 4. Numerical Differentiation and Integration

# Three-Point Midpoint Formula
# INPUT: f; x0; h
# OUTPUT: Approximation of Differentiation of f at x0 
def three_midpoint_differentiate(f, x0, h):
        return (f(x0+h)-f(x0-h))/(2*h)

# Three-Point Endpoint Formula
# INPUT: f; x0; h
# OUTPUT: Approximation of Differentiation of f at x0 
def three_endpoint_differentiate(f, x0, h):
        return (-3*f(x0)+4*f(x0+h)-f(x0+2*h))/(2*h)

# Five-Point Midpoint Formula
# INPUT: f; x0; h
# OUTPUT: Approximation of Differentiation of f at x0 
def five_midpoint_differentiate(f, x0, h):
        return (f(x0-2*h) -8*f(x0-h)+8*f(x0+h)-f(x0+2*h))/(12*h)

# Five-Point Endpoint Formula
# INPUT: f; x0; h
# OUTPUT: Approximation of Differentiation of f at x0 
def five_endpoint_differentiate(f, x0, h):
        return (-25*f(x0)+48*f(x0+h)-36*f(x0+2*h)+16*f(x0+3*h)-3*f(x0+4*h))/(12*h)

# Composite Simpson's Rule (Composite Numerical Integration)
# To approximate the integral I = integral from a to b f(x)dx:
# INPUT: Endpoints a, b; even positive integer n.
# OUTPUT: Approximation XI to I.
def composite_simpson_integral(f, a, b, n):
        
        if n%2 == 0:
                h = (b-a)/n
                XI_0 = f(a) + f(b)
                XI_1 = 0
                XI_2 = 0

                for i in range(n-1):
                        X = a + i*h

                        if i%2==0:
                                XI_2 = XI_2+f(X)
                        else:
                                XI_1 = XI_1+f(X)
                XI = h*(XI_0 + 2*XI_2 +4*XI_1)/3
                return XI
        else:
                return 'n should be even positive integer'
        

# Composite Trapezoidal Rule (Composite Numerical Integration)
# To approximate the integral I = integral from a to b f(x)dx:
# INPUT: Endpoints a, b; even positive integer n.
# OUTPUT: Approximation XI to I.
def composite_trapezoid_integral():
        return 1

def composite_midpoint_integral():
        return 1

# Romberg Integration
# To approximate the integral I = integral from a to b f(x)dx, select an integer n > 0.
# INPUT: Endpoints a, b; integer n.
# OUTPUT: An array R (Compare R by rows; only teh last two rows are saved in storage).
test_f = lambda x: x**2
def romberg_integration(f, a, b, n):
        
        h = b - a
        R = [ [ 0 for i in range(n) ] for j in range(n) ]
        
        R[0][0] = (h/2)*(f(a)+f(b))

        return R

def summation(f, x, n):
        temp_sum = 0        
        for i in range(n):
                temp_sum += f(x)
        return temp_sum

# Adaptive Quadrature
# To approximate the integral I = integral from a to b f(x)dx to within a given tolerance:
# INPUT: Endpoints a,b; tolerance TOL; limit N to number of levels.
# OUTPUT: Approximation APP or message that N is exceeded.
def adaptive_quadrature_integral():
        return 1

# Simpson's Double Integral
# To approximate the integral I = integral from a to b integral from c(x) to d(x) f(x,y) dy dx:
# INPUT: Endpoints a, b; even positive integers m,n.
# OUTPUT: Approximation J to I.
def simpson_double_integral(f, a, b, c, d, m, n):
        
        if m%2 == 0 and n%2 == 0:
                h = (b-a)/n
                J_1 = 0; J_2 = 0; J_3 = 0

                for i in range(n):
                        x = a + i*h
                        HX = (d-c)/m
                        K_1 = f(x,c)+f(x,d)
                        K_2 = 0; K_3 = 0

                        for j in range(1,m-1):
                                y = c + j*HX
                                Q = f(x,y)

                                if j%2==0:
                                        K_2 = K_2 + Q
                                else:
                                        K_3 = K_3 + Q
                                        
                        L =  ((K_1+2*K_2+4*K_3)*HX)/3
                        if i==0 or i==n:
                                J_1 = J_1 + L
                        elif i%2==0:
                                J_2 = J_2 + L
                        else:
                                J_3  = J_3 + L
                J = h*(J_1+2*J_2+4*J_3)/3
                return J
        else:
                return 'n should be even positive integer'

# Gaussian Double Integral
# To approximate the integral I = integral from a to b integral from c(x) to d(x) f(x,y) dy dx:
# INPUT: Endpoints a, b; even positive integers m,n.
# (The roots ri,j and coefficients ci,j need to be available for i = max{m,n} and for 1 <= j <= i)
# OUTPUT: Approximation J to I.
def gaussian_double_integral():
        return 1

# Gaussian Triple Integral
# To approximate
# the integral I = integral from a to b integral from c(x) to d(x) integral from e(x) to f(x) f(x,y,z) dz dy dx:
# INPUT: Endpoints a, b; even positive integers m,n.
# (The roots ri,j and coefficients ci,j need to be available for i = max{m,n} and for 1 <= j <= i)
# OUTPUT: Approximation J to I.
def gaussian_triple_integral():
        return 1

### 5. Initial-Value Problems for Ordinary Differential Equations

## Euler's Method
# y(t_i+1) = y(t_i) +h*f(t_i, y(t_i)).
# w0 = y0; w_i+1 = w_i + h*w(t_i, w_i), for each i = 0,1,...,N-1.
# To approximate the solution of the initial-value problem, dy/dt = f(t,y), a <= t <= b, y(a)=w0
# at (N+1) equally spaced numbers in the interval [a,b]:
# INPUT: Differential equation f(t,y); endpoints a, b; integer N; initial condition y0.
# OUTPUT: Approximation w to y at the (N+1) values of t.
def euler_method(f, a, b, N, y0):

    h = (b-a)/N
    t = a
    w = y0

    for i in range(1,N+1):
        w = w+h*f(t,w)
        t = a + i*h
        print('t_'+str(i)+': ', t, 'w_'+str(i)+': ', w)
        
    return w

# Runge-Kutta (Order Four)
# To approximate the solution of the initial-value problem y'=f(t,y), a <= t <= b, y(a) = ⍺,
# at (N+1) equally spaced numbers in the interval [a,b]:
# INPUT: Endpoints a, b; integer N; initial condition ⍺.
# OUTPUT: Approximation w to y at the (N+1) values of t.
f = lambda t, y : y-t**2+1
def runge_kutta(f,a,b,N,y0):
        h = (b-a)/N
        t = a; w = y0
        print('Initial Value (t0,y0) = ',t,',',w)

        for i in range(1, N):
                K_1 = h*f(t,w)
                K_2 = h*f(t+h/2, w+K_1/2)
                K_3 = h*f(t+h/2, w+K_2/2)
                K_4 = h*f(t+h, w+K_3)

                t = a + i*h
                w = w+(K_1+2*K_2+2*K_3+K_4)/6
        
        return t, w

# Runge-Kutta-Fehlberg Method
# To approximate the solution of the initial-value problem y'=f(t,y), a <= t <= b, y(a) = ⍺,
# with local truncation error within a given tolerance:
# INPUT: Endpoints a, b; integer N; initial condition ⍺; tolerance TOL; maximum step size hmax;
# minimun step size hmin.
# OUTPUT: t, w, h where w approximates y(t) and the step size h was used or a message that
# teh minimum step size was exceeded.
def runge_kutta_fehlberg():
        return 1

### 6. Direct Methods for Solving Linear Systems









    
