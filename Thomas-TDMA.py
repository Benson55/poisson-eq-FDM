#Solution of the Poisson equation within the domain [0,1] using
#Discretization is done both using second order Finite Difference schemes
#Dirichlet boundary conditions is implemented on both the boundaries
#The resulting linear algebraic system is solved using TDMA (Thomas algorithm)

#This is a part of assignments for the course "Computational Heat & Fluid Flow (ME 605)" taught at IIT Goa during the winter semester of 2020

#Author: Nithin Adidela, IIT Goa
#email: nithin.adidela.16003@iitgoa.ac.in
#Tested with python 3.7.5 on 12 February 2020

# ========== Importing Libraries ==========

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

# ********** Thomas algorithm to solve a tridiagonal matrix system **********
# Input diagonal elements are loaded after applying the bounday conditions (atleast in this code)
# [Input ] n -- number of mesh points
# [Input ] b -- main diagonal elements of the matrix
# [Input ] c -- lower diagonal elements of the matrix
# [Input ] a -- upper diagonal elements of the matrix
# [Input ] d -- Right Hand Side of the matrix system
# [Input ] x_thomas -- solution output
def Thomas(n, b, c, a, d, x_thomas):
    bdash = np.zeros(n)
    bdash[0]=b[0]/b[0]
    cdash=np.zeros(n)
    cdash[0]=c[0]
    ddash=np.zeros(n)
    ddash[0]=d[0]/b[0]
    for i in range(1, len(b)-1):
        bdash[i]=1
        cdash[i]=c[i]/(b[i]-a[i]*cdash[i-1])
        ddash[i]=(d[i]-a[i]*ddash[i-1])/(b[i]-a[i]*cdash[i-1])
    bdash[len(b)-1]=1
    ddash[len(b)-1]=(d[len(b)-1]-a[len(b)-1]*ddash[len(b)-2])/(b[len(b)-1]-a[len(b)-1]*cdash[len(b)-2])
    x_thomas[n-1]=ddash[n-1]/bdash[n-1]
    for i in range(n-2,-1,-1):
        x_thomas[i]=(ddash[i]-cdash[i]*x_thomas[i+1])/bdash[i]


# ========== Program begins here ==========
# ---------- Set the following input parameters  ----------
l=1.0			# length of the domain
n=11			# number of mesh divisions (for first question)	
del_x=l/(n-1)		# mesh size
x=np.zeros(n)		# position vector
bc       = 'Dirichlet'	# Specify the boundary conditions

# ---------- Initialize empty arrays  ----------
exact_sol=np.zeros(n)	# exact solution
c=np.zeros(n)		# lower diagonal
b=np.zeros(n)		# main diagonal
a=np.zeros(n)		# upper diagonal
d=np.zeros(n)		# Right Hand Side
sol_thomas=np.zeros(n)		# Solution from Thomas Algorithm
pi=np.pi		# A variable pi with the value Pi(π) using NumPy

# ---------- One for loop to append values to various arrays initialised above ----------
for i in range(1, len(x)-1): # for all the elements except the extremes
    x[i]=i*del_x	# to position vector
    a[i]=-1.0		# to upper diagonal        
    b[i]=2.0		# to main diagonal
    c[i]=-1.0		# to lower diagonal
    d[i]=(pi**2)*np.sin(pi*x[i])*(del_x**2)	# to Right Hand Side
    exact_sol[i]=np.sin(pi*x[i])		# to exact solution  	

x[0]=0		# left most value of position vector
x[len(x)-1]=1	# right most value of position vector

exact_sol[0]=np.sin(pi*x[0])			# left most value of exact solution
exact_sol[len(x)-1]=np.sin(pi*x[len(x)-1])	# right most value of exact solution

# ---------- Applying boundary conditions by changing the tridiagonal matrix extreme elements ----------
# ---------- Left hand boundary (Dirichlet) ----------
a[0]=0.0        
b[0]=1.0
c[0]=0.0
d[0]=0.0
# ---------- Right hand boundary (Dirichlet)----------
a[n-1]=0.0        
b[n-1]=1.0
c[n-1]=0.0
d[n-1]=0.0


# ---------- Calling the Thomas algorithm function to solve the tridiagonal matrix system ----------
error_Thomas = 0.0	# Initialise a variable for error produced whilw solving using Thomas algorithm	
Thomas(n, b, c, a, d, sol_thomas) # Calling the function
print(sol_thomas)		# Printing the solution array
error_Thomas = abs(max(sol_thomas-exact_sol, key=abs))
# Calculationg the error value as maximum of the absolute difference of elements in numerical, exact solutions	
print(error_Thomas)	# Printing the error

# ---------- Commands to plot the desired results ----------

plt.plot(x,exact_sol, label='Exact solution', color='green', marker='o', markerfacecolor='green', linestyle='dashed')
plt.plot(x,sol_thomas, label='Thomas Algorithm', color='blue', marker='x', markerfacecolor='blue', linestyle='dashed')
plt.title('Numerical Solution vs Exact Solution for 11 grid points')
plt.xlabel('Domain')
plt.ylabel('Numerical Solution vs Exact Solution')
plt.legend()  
plt.show() 

# ========== End of program ==========
