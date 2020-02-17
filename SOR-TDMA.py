#Solution of the Poisson equation within the domain [0,1] using
#Discretization is done both using second order as well as fourth order Finite Difference schemes
#Dirichlet boundary conditions is implemented on both the boundaries
#The resulting linear algebraic system is solved using Successive Over-relaxation algorithm with Over-relaxation parameter 1.5

#This is a part of assignments for the course "Computational Heat & Fluid Flow (ME 605)" taught at IIT Goa during the winter semester of 2020

#Author: Nithin Adidela, IIT Goa
#email: nithin.adidela.16003@iitgoa.ac.in
#Tested with python 3.7.5 on 12 February 2020

# ========== Importing Libraries ==========

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

# ********** This function produces a tridiagonal matrix **********
# In this code this function is called when Gauss-Jacobi, Gauss-Seidel and SOR methods are used but not for Thomas
# Input diagonal elements are loaded after applying the bounday conditions
# [Input ] l -- lower diagonal elements of the matrix
# [Input ] m -- main diagonal elements of the matrix
# [Input ] u -- upper diagonal elements of the matrix
def tridiag(l, m, u, k1=-1, k2=0, k3=1):
    l=np.zeros(n-1)
    for i in range(0, len(l)-1):
        l[i]=c[i+1]
    u=np.zeros(n-1)
    for i in range(0, len(u)):
        u[i]=a[i]
    return np.diag(l, k1) + np.diag(m, k2) + np.diag(u, k3)

# ********** Successive Over-Relaxation algorithm to solve a tridiagonal matrix system ***********
# Input diagonal elements are loaded after applying the bounday conditions (atleast in this code)
# [Input ] n -- number of mesh points
# [Input ] A_jac -- tridiagonal matrix
# [Input ] sol_jac -- solution matrix
# [Input ] B_jac -- Right Hand Side of the matrix system
# [Input ] iter_jac -- iteration count
# [Input ] norm_jac -- norm value
def SOR_TDM(n, A_SOR, X_SOR , B_SOR , iter_SOR, norm_SOR, w):
    x_SOR=np.zeros(n)
    Ax_SOR=np.zeros(n)
    converged = False
    while converged==False:
          for i in range(0, len(x_SOR)-1):
              if i==0:
                 x_SOR[i]= (1-w)*X_SOR[i]+w*(B_SOR[i]-A_SOR[[i],[i+1]]*X_SOR[i+1])/A_SOR[[i],[i]]
              elif i==n-1:
                 x_SOR[i]= (1-w)*X_SOR[i]+w*(B_SOR[i]-A_SOR[[i],[i-1]]*x_SOR[i-1])/A_SOR[[i],[i]]
              else: 
                 x_SOR[i]= (1-w)*X_SOR[i]+w*(B_SOR[i]-A_SOR[[i],[i-1]]*x_SOR[i-1]-A_SOR[[i],[i+1]]*X_SOR[i+1])/A_SOR[[i],[i]]
          Ax_SOR = np.dot(A_SOR,x_SOR) 
          norm_SOR = abs(max(B_SOR-Ax_SOR, key=abs))
          for i in range(0, len(x_SOR)):
              X_SOR[i]=x_SOR[i]
          iter_SOR = iter_SOR + 1
          print(iter_SOR, norm_SOR)
          if norm_SOR <0.00000001:
             converged = True 

# ********** This function applies boundary conditions on the given array **********
def apply_bc(u):
    u[0] = 0
    if bc == 'Dirichlet':
       u[n-1] = 0 
    else:
        sys.exit('Set correct boundary condition at the right side boundary')

# ========== Program begins here ==========
# ---------- Set the following input parameters  ----------
l=1.0			# length of the domain
n=11			# number of mesh divisions (for first question)	
del_x=l/(n-1)		# mesh size
x=np.zeros(n)		# position vector
bc       = 'Dirichlet'	# Specify the boundary conditions
L_inf_SOR=0		# L-infinity Norm (for second question)

# ---------- Initialize empty arrays  ----------
exact_sol=np.zeros(n)	# exact solution
c=np.zeros(n)		# lower diagonal
b=np.zeros(n)		# main diagonal
a=np.zeros(n)		# upper diagonal
d=np.zeros(n)		# Right Hand Side
sol_SOR=np.zeros(n)	# Solution from SOR
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

# ---------- Calling the function which produces a nxn tridiagonal matrix ----------
TDM = tridiag(c, b, a)

# ---------- Calling the Successive Over-relaxation algorithm to solve the tridiagonal matrix system ----------
norm_sor=1	
iter_sor=0
w=1.5		# Over-relaxation parameter		
SOR_TDM(n, TDM, sol_SOR , d , iter_sor, norm_sor, w)	# Successive Over-Relaxation
#print(norm_sor, iter_sor)	# Print norm as well as iteration count
#L_inf_SOR= abs(max(sol_SOR-exact_sol, key=abs))
#print(L_inf_SOR)

# ---------- Commands to plot the desired results ----------

plt.plot(x,exact_sol, label='Exact solution', color='green', marker='o', markerfacecolor='green', linestyle='dashed')
plt.plot(x,sol_SOR, label='Successiver Over-Relaxation with w=1.5', color='black', marker='s', markerfacecolor='black', linestyle='dashed')
plt.title('Numerical Solution vs Exact Solution for 11 grid points')
plt.xlabel('Domain')
plt.ylabel('Numerical Solution vs Exact Solution')
plt.legend()  
plt.show() 

# ========== End of program ==========
