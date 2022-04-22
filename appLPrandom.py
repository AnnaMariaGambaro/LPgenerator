# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:10:08 2021

@author: Anna
"""
import streamlit as st
import LPrandom as myfn

#option = st.sidebar.selectbox("What do you want?", ('Create new LP problem'))

st.write(""" # Linear Program Generator """) 
st.write('Create random feasible linear programs and solve them')     
st.latex(r'''
\begin{aligned}
    \max & \: c^T x \\
    \textrm{sub to } & \: A x \leq b \\
    &  x \geq 0 \\
\end{aligned}
''')


# An alias for our state
state = st.session_state

# A function to easily go from one step to another
def change_step(next_step):
    state.step = next_step

# Let's initialize our session state
if "data" not in state:
    state.data = []
    state.step = "init"

# Step 1
if state.step == "init":
    st.button("Create new LP", on_click=change_step, args=["create"])

# Step 2
if state.step == "create":
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

    A = LPproblem['A']
    b = LPproblem['b']   

    col1, col2, col3 = st.columns(3)       
    with col1:                          
        st.write('Constraint matrix A = ', A) 
    with col2:
        st.write('Constraint vector b  = ', b)
    with col3:
        st.write('Obj. fun. coeff. c = ', c) 
    
    state.data = LPproblem
    st.button("Show solution", on_click=change_step, args=["solution"])


if state.step == "solution": 
    st.button("Create new LP", on_click=change_step, args=["create"])
    
    LPproblem = state.data
    A = LPproblem['A']
    b = LPproblem['b'] 
    c = LPproblem['c']          
    col1, col2, col3 = st.columns(3)       
    with col1:                          
        st.write('Constraint matrix A = ', A) 
    with col2:
        st.write('Constraint vector b  = ', b)
    with col3:
        st.write('Obj. fun. coeff. c = ', c)     
    
    sol = LPproblem['sol'] 
    OF = LPproblem['obj_fun']
    vertex = myfn.vertex(LPproblem)
    
    col1, col2, col3 = st.columns(3)       
    with col1:                          
        st.write('Solution x = ', sol) 
    with col2:
        st.write('Obj. Fun. = ', OF)
    with col3:
        st.write('Vertices = ', vertex) 
    
    fig = myfn.regiongraph(LPproblem, 500)
    fig
    