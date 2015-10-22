# Numerical Method (Analysis) Library - pynumerical.py
# Date: Oct/21/2015
# Author: Minwoo Bae
import math

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
        	
def bisection_method(f_x, a_n, b_n, num):

        def isNegative(f_x, a_n, b_n):
                if f_x(a_n) * f_x(b_n) < 0: 
                        return True
                else: 
                        return False
        
        def findPivot(a_n, b_n):
                return (a_n+b_n)/2
        
	i = 0
	pivot = 0
	if isNegative(a_n, b_n):
		while i < num:
			pivot = findPivot(a_n, b_n)
			if f_x(pivot) > 0 and f_x(b_n) < 0:
				a_n = pivot
				print('Pivot (=>a) ['+str(i)+'] = '+str(pivot))
				print('F(Pn) = '+str(f_x(pivot))) 
				print('')
			else:
				b_n = pivot
				print('Pivot (=>b) ['+str(i)+'] = '+str(pivot))
				print('F(Pn) = '+str(f_x(pivot)))
				print('')
			i += 1
	else: 
		print('Please reset up the right interval.')
	return pivot

def nevilles_method(x0, x, y):

	n = len(x)

	# Create and initiate an n x n multidimentional array: 
	q = [ [ 0 for i in range(n) ] for j in range(n) ]

	# Set y[i] value into q[i][0]:
	for i in range(n):
		q[i][0] = y[i]
	
	for i in range(1, n):
		for j in range(1, i):
			q[i][j] = ((x0 - x[i - j])*(q[i][j - 1]) - (x0 - x[i])*(q[i - 1][j - 1]))/(x[i] - x[i - j])

	print("Result of Neville's Method:\n");
	for i in range(n):
		for j in range(n):
			print("%.4f" %q[i][j])
		print("\n")
