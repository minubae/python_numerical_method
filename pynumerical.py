# Numerical Method (Analysis) Module - pynumerical.py
# Date: Oct/21/2015, Wed - 
# Version: ver 1.0
# Author: Minwoo Bae
# Contact: minubae.nyc@gmail.com

import math

#######################################################################################
########################  A Function for the purpose of testing #########################

fx = lambda x: 2**x - 3
fx2 = lambda x:  x**3 + 4*x**2 - 10
fx3 = lambda t,y: y/t - (y/t)**2  # --> Differential Equation; for testing of this Euler's method.

#########################################################################################

#########################################################################################
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

################################### Fixed-Point Iteration ##################################
# Fixed-Point Theorem:
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

# def fixed_point(fx, p, num):
#     i = 1
#     p_n = p
#     fx(p_n)
#     while i < num: 
#         print('P_'+str(i)+': '+str(fx(p_n)))
#         p_n = fx(p_n)
#         i+=1
#     return p_n


################################### Newton's Method ##################################
def newthon_method(p_0, tol, num):
	i = 0
	while i < num:
		p = p_0 - f(p_0)/fPrime(p_0)
		if math.fabs(p-p_0) < tol:
			return [p, i]; break
		i+=1
		p_0 = p

	print('The method failed after N_0 interations, N_0 = '+str(num))

################################### Secant Method ##################################
def secant_method(p_0, p_1, tol, num):
	i = 2
	while i <= num:
		p = p_1 - (f(p_1)*(p_1 - p_0))/((f(p_1) - f(p_0)))
		if math.fabs(p - p_1) < tol:
			return [p, i]
			break
		i+=1
		p_0 = p_1
		p_1 = p

	print('The method failed after N_0 iterations, N_0 = '+str(num))

################################### The Method of False Postion ##################################
def false_position(p_0, p_1, tol, num):
	i = 2
	q_0 = f(p_0)
	q_1 = f(p_1)
	while i <= num:
		p = p_1 - q_1*(p_1-p_0)/(q_1 - q_0)
		if math.fabs(p - p_1) < tol:
			return [p, i]
			break
		i+=1
		q = f(p)
		if q*q_1 < 0:
			p_0 = p_1
			q_0 = q_1
		p_1 = p
		q_1 = q
	print('The method failed after N_0 iterations, N_0 = '+str(num))
	


##################################### Neville's Method ######################################
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

############################ Newton's Divided-Difference Formula ############################
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

############################ Hermote Interpolation ##############################
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



################################  Euler's Method ###################################
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











    
