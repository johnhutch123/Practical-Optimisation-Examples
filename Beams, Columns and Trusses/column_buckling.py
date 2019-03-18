# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 08:42:02 2019

Column bluckling optimisation using analytical equations.

See page 10 of 'Engineering Optimization: Theory and Practice' by Singiresu S. Rao

@author: John Hutcheson
"""

from scipy import optimize
import numpy as np
from matplotlib import pyplot as plt


# Define equations for objective and constraints.
def obj(X):
    return (5*rho*(np.pi*X[0]*X[1]*L))+(2*X[0])

def sig_1(X):
    return sig_yield - F/(np.pi*X[0]*X[1])

def sig_2(X):
    return ((np.pi**2)*E/(8*L**2)*(X[0]**2+X[1]**2)) - F/(np.pi*X[0]*X[1])

def callback(X):
    history.append(X)

if __name__ == '__main__':

    
    # Define the problem
    F = 24500
    E = 8.34e+10
    rho = 2500
    L = 2.5


    # Constraint limits
    sig_yield = 4.91e7
    
    dia_min = 0.02
    dia_max = 0.14
    
    th_min = 0.002
    th_max = 0.008
    
    bound_min = np.array([dia_min, th_min])
    bound_max = np.array([dia_max, th_max])
    
    
    # Set up the optimisation
    X_init = np.array([dia_max,th_max]) # Starting point at largest diameter and thickness.
    
    # Save the design variables.
    history = [X_init] # initalise with intial starting point.
    
    constraints = [{'type':'ineq','fun':sig_1},
                   {'type':'ineq','fun':sig_2}]
     
    bounds = optimize.Bounds(bound_min,bound_max)
    
    
    # Run the optimisation
    result = optimize.minimize(obj,X_init,method='SLSQP',bounds=(bounds),constraints=(constraints),callback=(callback))
    print(result)
    
    
    # Visualise results 
    history = np.array(history)
    
    plt.plot(history[:,0],history[:,1],'.-')
    plt.xlabel('Diameter, (m)')
    plt.ylabel('Thickness, (m)')
    plt.title('Convergence History')
    plt.grid(True)
    
    thConst = np.linspace(th_min,th_max,20)
    
    # yield stress constraint
    dConst1 = F / (sig_yield * np.pi * thConst)
    
    
    # buckling stress constraint
    
    
    # Plot constraint curves
    
    plt.fill_between(dConst1,thConst,color='r')
    
    # Yield stress constraint
    
    
    # Buckling stress constraint
