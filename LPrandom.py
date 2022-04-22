# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:15:23 2021

@author: Anna
"""
import numpy as np
import random
from scipy.optimize import linprog
import matplotlib.pyplot as plt

"""
m = 3
n=2
muA = 0 
sigmaA = 1
rho = 0.8 
mu_alpha = 0 
sigma_alpha = 1 
gamma = 0.7
"""

def matrixgen(m, n, muA, sigmaA, rho):
    tol = 0.001
    A = np.zeros((m,n))
    segno = np.random.randint(2, size=(m,n))*2-1
    if type(muA) == int or type(muA) == float:
        temp = np.linspace(0,m-1,m)
        muA = temp*0+muA        
    for ii in range(m):
        for jj in range(n):
            r = random.random()
            if r<rho:
                A[ii,jj] = np.random.normal(muA[ii], sigmaA, 1)
            else:
                A[ii,jj]= 0 
    if np.abs(A).max() < tol:
        for jj in range(n):
            A[0,jj] = np.random.normal(muA[ii], sigmaA, 1)
    for ii in range(m):
        if np.abs(A[ii,:]).max() < tol:
            for jj in range(n):
                A[ii,jj] = np.random.normal(muA[ii], sigmaA, 1)    
    A = segno*A
    return A


def encodegen(m, n, muA, sigmaA, rho, mu_alpha, sigma_alpha, gamma ):
    A = matrixgen(m, n, muA, sigmaA, rho)
    alpha = np.random.lognormal(mu_alpha, sigma_alpha, m+n)
    k = round(gamma*min(m,n))
    k2 = round(gamma*n)
    beta = np.zeros((1,m+n))
    indices = np.random.choice(n, k, replace = False)
    indices2 = np.random.choice(m, m-k, replace = False)
    indices3 = np.random.choice(n, k2, replace = False)
    for ii in range(n):
        if ii in indices:
            #print(ii)
            beta[0,ii] = 1
    for jj in range(m):
        if jj in indices2:
            #print(jj+n)
            beta[0,jj+n] = 1
    beta = beta[0,:]
    """
    for ii in range(n):
        if ii in indices3:
            aux = np.random.beta(0.5,0.5)
            alpha[ii] = max(alpha[ii].round(0),1)-aux
        else:
            alpha[ii] = max(alpha[ii].round(0),1)
    """
    return A, alpha, beta

def LPconstructor(m,n, muA =None, sigmaA = None, rho = None, mu_alpha = None, 
                  sigma_alpha= None, gamma= None, digit = None):
    if muA == None:
        muA = 0
    if sigmaA == None:
        sigmaA = 1
    if rho == None:
        rho = 0.8
    if mu_alpha == None:
        mu_alpha = 0
    if sigma_alpha == None:
        sigma_alpha = 1
    if gamma == None:
        gamma = 0.8
    if digit == None:
        digit = 2
    
    encode = encodegen(m, n, muA, sigmaA, rho, mu_alpha, sigma_alpha, gamma )
    A = encode[0].round(digit)
    alpha = encode[1].round(digit)
    beta = encode[2]
    
    sol =  beta[0:n]*alpha[0:n]
    slack_dual = (1-beta[0:n])*alpha[0:n]
    sol_dual = (1-beta[n:m+n])*alpha[n:m+n]
    slack = beta[n:m+n]*alpha[n:m+n]
    b = np.dot(A,sol)+slack
    c = np.dot(A.transpose(),sol_dual)-slack_dual
    
    fo = np.dot(c,sol)
    
    res = dict()
    res['A'] = A
    res['b'] = b
    res['c'] = c
    res['sol'] = sol
    res['sol_dual'] = sol_dual
    res['slack'] = slack
    res['slack_dual'] = slack_dual
    res['obj_fun'] = fo
    return res

def checkLPsol(m,n, muA =None, sigmaA = None, rho = None, mu_alpha = None, 
               sigma_alpha= None, gamma= None, digit = None):
    if muA == None:
        muA = 0
    if sigmaA == None:
        sigmaA = 1
    if rho == None:
        rho = 0.8
    if mu_alpha == None:
        mu_alpha = 0
    if sigma_alpha == None:
        sigma_alpha = 1
    if gamma == None:
        gamma = 0.8
    if digit == None:
        digit = 2
    
    res = LPconstructor(m,n, muA, sigmaA, rho, mu_alpha, sigma_alpha, gamma, digit)
    A = res['A']
    b = res['b']
    c = res['c']
    
    #-c insted of c because linprog minimize instead of maximize
    sol = linprog(-c, A_ub=A, b_ub=b)
    tol = 0.000001 # tollerance
    if abs(res['obj_fun']+sol.fun) < tol:
        check = True
    else:      
        check = False      
        
    return check, res, sol

def checkpoint(LPproblem, x):
    A = LPproblem['A']
    b = LPproblem['b']
    dim = A.shape
    temp=1
    tol = 0.0000001
    for kk in range(dim[0]):
        temp = temp*(np.dot(A[kk,:],x) <= b[kk]+tol )
    temp = temp*(x.min()>=0) 
    return temp 

def vertex(LPproblem):
    A = LPproblem['A']
    b = LPproblem['b']
    sol = LPproblem['sol']
    
    dim = A.shape
    if dim[1] > 2:
        print('Errore: più di due variabili')
        return 0
    else:
        vertex = np.copy(sol)
        for ii in range(dim[0]):
            if abs(A[ii,1]) > 0:
                x = np.array([0, b[ii]/A[ii,1]])
                condition = checkpoint(LPproblem, x)
                if condition == True:
                    vertex = np.concatenate((vertex, x))
            if abs(A[ii,0]) > 0:
                x = np.array([b[ii]/A[ii,0],0])  
                condition = checkpoint(LPproblem, x)
                if condition == True:
                    vertex = np.concatenate((vertex, x))
            for jj in range(ii+1,dim[0]):
                try:
                    mat = A[[ii,jj],:]
                    vet = b[[ii,jj]]
                    x = np.linalg.solve(mat,vet)
                    condition = checkpoint(LPproblem, x)
                    if condition == True:
                        vertex = np.concatenate((vertex, x))
                except:
                    vertex = vertex
                    #print('warning: vincolo (%s' %ii + ',%s)' %jj)
        vertex = vertex.reshape(round(vertex.shape[0]/2),2)
        vertex = np.unique(vertex.round(3),axis=0)
        return vertex


def regiongraph(LPproblem,nstep = None):
    if nstep == None:
        nstep = 100
    
    vert = vertex(LPproblem)
    A = LPproblem['A']
    b = LPproblem['b']
    c = LPproblem['c']
    sol = LPproblem['sol']
    
    dim = A.shape
    if dim[1] > 2:
        print('Errore: più di due variabili')
    else:
        d_max = vert.max()+1
        d_min = max(0,vert.min()-1)
        """
        #find maximum values for horizontal axis
        if abs(A[:,0]).min() == 0 and abs(A[:,1]).min() == 0:
            d_max = max(sol[0], sol[1])
        elif abs(A[:,1]).min() == 0:  
            aux = abs(b/A[:,0])
            d_max = max(max(sol[0], sol[1]),aux.max())
        elif abs(A[:,0]).min() == 0:
            aux = abs(b/A[:,1])
            d_max = max(max(sol[0], sol[1]),aux.max())   
        else:
            aux = abs(b/A[:,0])
            aux2 = abs(b/A[:,1])
            d_max = max(max(sol[0], sol[1]),aux.max(), aux2.max())     
            
        if d_max == 0:
            d_max = A.max()
        """
        
        fig = plt.figure(figsize=(9,9))
        # plot the constraints
        d = np.linspace(d_min,d_max,nstep)
        x = d
        y = x*0
        plt.plot(x, y,'k')
        y = x
        x = y*0
        plt.plot(x, y,'k')
        for jj in range(dim[0]):
            if abs(A[jj,1]) > 0:
                x = d
                y = b[jj]/A[jj,1] - (A[jj,0]/A[jj,1])*x
                leg = '%s x + ' %A[jj,0] + '%s y' %A[jj,1] + ' <= %s' %b[jj]
                plt.plot(x, y, label = leg)
            else:
                y = d
                x = y*0+b[jj]/A[jj,0]
                leg = '%s x + ' %A[jj,0] + '%s y' %A[jj,1] + ' <= %s' %b[jj]
                plt.plot(x, y, label = leg)
        # plots the feasible region
        x,y = np.meshgrid(d,d)
        temp = x*0+True
        for jj in range(dim[0]):
            temp = np.logical_and(temp, A[jj,0]*x+A[jj,1]*y<=b[jj])
        plt.imshow(temp.astype(int), extent=(x.min(),x.max(),y.min(),y.max()), 
                   origin="lower", cmap="Greys", alpha = 0.3)
        
        #plot the objective function
        fopt = np.dot(c,sol)
        #print('FO = %s' %fopt)
        x = d
        y = (-c[0]/c[1])*x + fopt/c[1]
        plt.plot(x, y, '--',label = 'F. O. = %s x ' %c[0] + '+ %s y' %c[1])
        plt.legend()
        plt.xlim([0,vert[:,0].max()+1])
        plt.ylim([0,vert[:,1].max()+1])
        return fig       
        
 

    
