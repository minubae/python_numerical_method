# Numerical Method (Analysis) Module - pynumerical.py
# Date: Oct/21/2015, Wed - 
# Version: ver 1.0
# Author: Minwoo Bae
# Contact: minubae.nyc@gmail.com

import math

#########################################################################################
########################  A Function for the purpose of testing #########################
def fx(x):
	return 2**x - 3
#########################################################################################

#########################################################################################
# Suppose that a_p is an approximation to p.
# The actual error is p-a_p,
# the abosulte error is |p-a_p|,
# and the relative error is |p-a_p|/|p|,provided that p != 0.
def absolute_value(p):
        return math.fabs(p)

def actual_error(p,a_p):
        return p - a_p

def absolute_error(p,a_p):
        return math.fabs(p-a_p)

def relative_error(p,a_p):
        result = 0
        if absolute_value(p) != 0:
                result = absolute_error(p,a_p)/absolute_value(p)
                return result
        else:
                print('error: absolute_value must not be equal to ',absolute_value(p))
##########################################################################################

############################ Bisection Method (Binary Search) ############################
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
def bisection_method(fx, a_n, b_n, num):

        i = 0
        pivot = 0
        
        def isNegative(fx, a_n, b_n):
                if fx(a_n) * fx(b_n) < 0: 
                        return True
                else: 
                        return False
        
        def findPivot(a_n, b_n):
                return (a_n+b_n)/2
        
        if isNegative(fx, a_n, b_n):
                for i in range(num):
                        pivot = findPivot(a_n, b_n)
                        if fx(pivot) > 0 and fx(b_n) < 0:
                                a_n = pivot
                        else:
                                b_n = pivot
                return pivot
        else:
                print('Please reset up the right interval.')

############################ Fixed-Point Iteration ############################               
def fixed_point(fx, p, num):
	i = 1
	p_n = p
	fx(p_n)
	while i < num: 
		print('P_'+str(i)+': '+str(fx(p_n)))
		p_n = fx(p_n)
		i+=1
	return p_n

############################ Neville's Method ############################
def neville_method(): # x0, y, fx
        x0 = 1.5
        x = [1.0, 1.3, 1.6, 1.9, 2.2]
        fx = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
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

############################ Newton's Divided-Difference Formula ############################
def divided_differences(): # x, fx
        x = [1.0, 1.3, 1.6, 1.9, 2.2]
        fx = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
        n = len(x)
        F = [ [ 0 for i in range(n) ] for j in range(n) ]
        
        for i in range(n):
                F[i][0] = fx[i]
        
        for i in range(1, n):
                for j in range(1, i+1):
                        F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (x[i] - x[i-j])
        return F

############################ Hermite Interpolation ############################