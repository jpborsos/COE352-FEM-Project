# John Paul Borsos
# COE 352: 1D Finite Element Project
# Backward Euler

import matplotlib.pyplot as plt 
import numpy as np
from numexpr import evaluate 
from math import pi, sqrt, e
from scipy.linalg import inv
global_expr = {'pi': pi, 'e': e}

### INPUTS / TERMS TO MODIFY: ###

# analytic solution
u_real = 'e**(-t)*sin(pi*x)'

# f(x,t)
f = '(pi**2-1)*e**(-t)*sin(pi*x)'

# mesh
x = np.linspace(0,1,11) # bottom of domain, top of domain, number of nodes
h = np.diff(x) # stepsize
numN = len(x) # number of nodes
numE = numN - 1 # number of elements

# dirichlet boundary conditions
u_0 = 'sin(pi*x)' # u(x,0)
dir_bc = np.array(([0, 0], [numN-1,0])) # in format [node, value]

# time
dt = 1/5 #timestep
time = np.arange(0,1+dt,dt) # total

### END OF INPUTS ###



# initialize matrices
M = np.zeros((numN,numN))    # global mass matrix
K = np.zeros((numN,numN))    # global stiffness matrix
K_loc = np.zeros((2,2))    # local stiffness matrix
M_loc = np.zeros((2,2))    # local mass matrix

# set gaussian quadrature values
gaussW = np.array([1., 1.])  # quadrature weights
gaussP = np.array([-1/sqrt(3), 1/sqrt(3)])  # quadrature points

# basis functions in eta space 
psi_1 = '(1-xi)/2'
psi_2 = '(1+xi)/2'
dpsi_1 = -1/2
dpsi_2 = 1/2
dpsi = np.array([-1/2, 1/2])

# evaluate basis functions at gaussian quadrature points
psi_1q = np.array(evaluate(psi_1, local_dict={'xi': gaussP}))
psi_2q = np.array(evaluate(psi_2, local_dict={'xi': gaussP}))
psi_q = np.vstack((psi_1q, psi_2q))

# create mapping
loc2glob_map = np.arange(numE)   
loc2glob_map = np.vstack((loc2glob_map, loc2glob_map+1)).T

# build M and K matrices
for k in range(numE):
    for l in range(2):
        for m in range(2):
            M_loc[l][m] = (h[k]/2) * np.sum(psi_q[l,:] * psi_q[m,:] * gaussW)
            K_loc[l][m] = 2* (2/h[k]) * dpsi[l] * dpsi[m]

    # combine local terms to get global matrix
    for l in range(2):
        global_node = loc2glob_map[k][l]
        for m in range(2):
            global_node2 = loc2glob_map[k][m]
            K[global_node][global_node2] += K_loc[l][m]
            M[global_node][global_node2] += M_loc[l][m]

# calculate inverses
B = (1/dt)*M + K
B_inv = inv(B)
B_inv_copy = np.copy(B_inv)
B_inv_M = B_inv @ M
B_inv_M_copy = np.copy(B_inv_M)

# initialize u by evaluating it at boundary condition
u = evaluate(u_0, local_dict={'x':x}, global_dict=global_expr)

for t in time:
    F = np.zeros(numN)

    for k in range(numE):
        x_quad = (h[k]/2) * (gaussP+1) + x[k] # map quadrature points back to x
        for l in range(2):
            F_quad = evaluate(f, local_dict={'t':t, 'x':x_quad}, global_dict=global_expr)
            F_loc = (h[k]/2) * np.sum(gaussW * F_quad * psi_q[l,:])
            global_node = loc2glob_map[k][l]
            F[global_node] += F_loc
    
    # impose boundary conditions
    for bc in dir_bc:
        j = bc[0]
        B_inv_M[j,:] = 0
        B_inv_M[:,j] = 0
        B_inv[j,j] = 1

    for j in range(numN):
        if j in dir_bc[:,0]:
            loc = np.where(dir_bc[:,0] == j)
            F[j] = dir_bc[loc,1]
            for i in range(numN):
                if i != j:
                    F[i] -= B_inv_copy[i,j]*dir_bc[loc,1]/dt 

    # perform backward euler to get next u vector
    u_new = (1/dt) * B_inv_M @ u + B_inv@F
    u = u_new


### PLOTTING
fig, ax = plt.subplots()

x_real = np.linspace(0,1,100)
u_real = evaluate(u_real, local_dict = {'x': x_real, 't': 1}, global_dict=global_expr)
ax.plot(x, u, label = 'Backward Euler')
ax.plot(x_real, u_real, label = 'Analytical')
ax.set_xlabel('x')
ax.set_ylabel('u')
ax.set_title('Backward Euler vs Analytical Solution @ t=1 (dt = 1/5)')
ax.legend()
ax.grid(True)

plt.savefig('./plots/Backward5.png',dpi = 200)
plt.show()