# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 08:33:40 2019

Truss optimisation using analytical equations.

See page 467 of 'Engineering Optimization: Theory and Practice' by Singiresu S. Rao

@author: John Hutcheson
"""

from scipy import optimize
import numpy as np
from matplotlib import pyplot as plt


# Define equations for objective and constraints.
def vol(A):
    return (20*(2**(1/2))*A[0])+(10*A[1])

def sig_1(A):
    return sig_up - (F*((2*(2**(1/2)))+(A[1]/A[0]))/(2*(2**(1/2))*A[0]+(1+2**(1/2))*A[1]))

def sig_2(A):
    return sig_up - F*(1/((A[0])+((2**(1/2))*A[1])))

def sig_3(A): 
    return (F*((2*2**(1/2))+(A[0]/A[1]))/(2*2**(1/2)*A[0]+(1+2**(1/2))*A[1])) - sig_low

def u_x(A):
    return disp_up - (F*L)/(E*A[0])*((2*A[0]+2**(1/2)*A[1]))/(2*2**(1/2)*A[0]+(1+2**(1/2))*A[1])

def u_y(A):
    return disp_up - (F*L)/(E*A[0])*(2*A[0])/(2*2**(1/2)*A[0]+(1+2**(1/2))*A[1])

def callback(X):
    history.append(X)

if __name__ == '__main__':

    
    # Define the problem
    F = 88964
    L = 0.254
    E = 6.89e+10


    # Constraint limits
    sig_up = 1.379e8
    sig_low = -1.034e8
    disp_up = 0.00508
    area_min = 0.000001
    area_max =  0.1
    area_start = 0.001
    
    
    # Set up the optimisation
    A_init = np.array([area_start,area_start]) # Starting point
    
    
    # Save progress of optimisation
    history = [A_init] # initalise with intial starting point.
    
    
    constraints = [{'type':'ineq','fun':sig_1},
                   {'type':'ineq','fun':sig_2},
                   {'type':'ineq','fun':sig_3},
                   {'type':'ineq','fun':u_x},
                   {'type':'ineq','fun':u_y}]
     
    bounds = optimize.Bounds(area_min,area_max)# place bounds on the memeber area to prevent it from going to zero.
    
    
    # Run the optimisation
    result = optimize.minimize(vol,A_init,method='SLSQP',bounds=(bounds),constraints=(constraints),callback=(callback))
    print(result)
    
    
    # Visualise the optimisation
    
    history = np.array(history)
    
    plt.plot(history[:,0],history[:,1],'.-')
    plt.xlabel('A1 = A3, (m^2)')
    plt.ylabel('A2, (m^2)')
    plt.title('Convergence History')
