# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:05:14 2021

@author: Anna
"""


import LPrandom as myfn

"""
print(myfn.checkLPsol(3,2))

print(myfn.checkLPsol(10,10))
"""

print('Cumstom parameters')
m = 3 #number of constraints
n=2 #number of variables, code run only for n=2
muA = [3,5,2] 
sigmaA = 2
rho = 0.6#0.99 # 1-rho = prob(A[ii,jj] = 0)
mu_alpha = 0 
sigma_alpha = 1 
gamma = 0.99
digit = 0

c_min = -1
OF = -1
dim = -1
ii = 0
while c_min <= 0 or OF <= 0 or dim < 0:
    res = myfn.checkLPsol(m,n, muA, sigmaA, rho, mu_alpha, sigma_alpha, gamma, digit)
    LPproblem = res[1] 
    c = LPproblem['c']
    OF = LPproblem['obj_fun']
    vertex = myfn.vertex(LPproblem)
    dim = vertex.shape[0]-2
    c_min = c.min()
    ii +=1
    #print(LPproblem['sol'])

print(res)
print('numero tentativi = %s' %ii)
"""
A = np.array([[ 11., -12.], [ -2.,  -4.], [ -2.,  -2.]])
b = np.array([-12.,  -3.,  -2.])
c = np.array([  9., -14.])
sol = np.array([0., 1.])


LPproblem = dict()
LPproblem['A'] = A
LPproblem['b'] = b
LPproblem['c'] = c
"""

grafico = myfn.regiongraph(LPproblem, 500)
print(vertex)