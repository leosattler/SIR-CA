###########################################################################
# SIR epidemiological model with bidimensional Celular Automata (CA).
# The reference for this implementation is:
# White, S. H., Del Rey, A. M., & Sánchez, G. R. (2007).
# Modeling epidemics using cellular automata. Applied
# mathematics and computation, 186(1), 193-202.
#
# INPE, Sao Jose dos Campos, SP, Brasil - July 6th, 2020
# Leonardo Sattler Cassara - leocassara@igeo.ufrj.br
###########################################################################

#==========================================================================
#                                 IMPORTS
#--------------------------------------------------------------------------
import numpy as np
from SIRca_output_handler import *
#==========================================================================
# END OF IMPORTS

#==========================================================================
#                                 INPUTS
#--------------------------------------------------------------------------
#                              Output folder
#--------------------------------------------------------------------------
# Folder must exist and be empty before run
output_folder = '.'
#--------------------------------------------------------------------------
#                            Grid parameters
#--------------------------------------------------------------------------
# Resolution of Grid (assuming a square grid)
resolution = 50
#
# Number of neighbors cells to consider during simulation.
# Muste be either 4 (Von Neumann) or 8 (Moore).
neighbors = 8
#
# Maximum iteration time
# Must be greater than 1
t_max = 150
#--------------------------------------------------------------------------
#                          Infection parameters
#--------------------------------------------------------------------------
# Population size (N)
#--------------------------------------------------------------------------
# Must be 1 or 2
case_N = 1
# where:
def set_N(grid, case_N):
    # Homegenous distribution of inhabitants
    # Setting to N_pop (population number per cell)
    if case_N == 1:
        N_pop = 100
        N[:,:] = N_pop
    # Increasing population from left to right
    if case_N == 2:
        for i in np.arange(l1, lpy2):
            for j in np.arange(c1, cpy2):
                grid[i, j] = int(np.exp(j))
    #
    return grid
#--------------------------------------------------------------------------
# Infected state
#--------------------------------------------------------------------------
# Must be 1
case_I = 1
# where:
def set_I(grid, case_I):
    # Central cell inital value only
    if case_I == 1:
        grid[int(np.shape(grid)[0]/2), \
             int(np.shape(grid)[1]/2)] = .3
    #
    return grid
#--------------------------------------------------------------------------
# Susceptible state
#--------------------------------------------------------------------------
# Must be 1
case_S = 1
# where:
def set_S(grid, case_S):
    # All individuals start with susceptibility = 1
    grid[:,:] = 1 
    # Central cell inital value only
    if case_S == 1:
        grid[int(np.shape(grid)[0]/2), \
             int(np.shape(grid)[1]/2)] = .7
    #
    return grid
#--------------------------------------------------------------------------
# Recovered state
#--------------------------------------------------------------------------
# Must be 1
case_R = 1
# where:
def set_R(grid, case_R):
    # Central cell inital value only
    if case_R == 1:
        grid[int(np.shape(grid)[0]/2), \
             int(np.shape(grid)[1]/2)] = 0.
    #
    return grid
#--------------------------------------------------------------------------
# Portion of infected individuals
# (that recover from the desease at each time step)
#--------------------------------------------------------------------------
# Must be between 0 and 1
e = .4
#--------------------------------------------------------------------------
# Connection factor
#--------------------------------------------------------------------------
# Must be 1 or 2
case_c = 1
# where:
def set_c(grid, case_c):
    # Constant connection factor
    if case_c == 1:
        # Setting to c_cte
        c_cte = 1.
        grid[:,:] = c_cte
    # Case of article example
    if case_c == 2:
        grid = connection_cond(grid)
    #
    return grid
#--------------------------------------------------------------------------
# Movement factor
#--------------------------------------------------------------------------
# Must be 1
case_m = 1
# where:
def set_m(grid, case_m):
    # Constant movment factor
    if case_m == 1:
        # Setting to c_cte
        m_cte = .5
        grid[:,:] = m_cte
    #
    return grid
#--------------------------------------------------------------------------
# Virulence of the epidemic
#--------------------------------------------------------------------------
# Must be between 0 and 1
v = .6
#--------------------------------------------------------------------------
# Vaccination rate
#--------------------------------------------------------------------------
# Must be between 0 and 1
omega = 0.
#==========================================================================
# END OF INPUTS

#==========================================================================
#                            CREATING GRIDS
#--------------------------------------------------------------------------
# Creating main grids 
#--------------------------------------------------------------------------
# Number of boundary ghost points (must be 1)
nghost = 1
# Population
N = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
# Infected
I = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
# Susceptible
S = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
# Recovered
R = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
# Grid of states
Q = np.array([I, S, R])
# Connection factor
c = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
# movement factor
m = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
#--------------------------------------------------------------------------
# Configuring horizontal ghost points
#--------------------------------------------------------------------------
l1 = nghost              # first physical poit
l2 = resolution+nghost-1 # last physical point
lpy2 = l2+1              # python's call for last physical
                         # (in fact first last ghost)
#--------------------------------------------------------------------------
# Configuring vertical ghost points
#--------------------------------------------------------------------------
c1 = nghost              # first physical poit
c2 = resolution+nghost-1 # last physical point
cpy2 = c2+1              # python's call for last physical
                         # (in fact first last ghost)
#==========================================================================
# END OF CREATING GRIDS

#==========================================================================
#                            AUXILIAR METHODS
#--------------------------------------------------------------------------
# Connection factor inital condition according to the article example
#--------------------------------------------------------------------------
def connection_cond(grid):
    half_grid = int(np.shape(grid)[0]/2)
    for i in np.arange(l1, lpy2):
        for j in np.arange(c1, cpy2):
            #
            if 0<=i and i<=half_grid and \
               0<=j and j<=half_grid:
                grid[i,j] = .6
            #
            if 0<=i and i<=half_grid and \
               (half_grid+1)<=j and j<=(2*half_grid-1):
                grid[i,j] = 1.
            #
            if (half_grid+1)<=i and i<=(2*half_grid-1) and \
               0<=j and j<=half_grid:
                grid[i,j] = 0.
            #
            if (half_grid+1)<=i and i<=(2*half_grid-1) and \
               (half_grid+1)<=j and (2*half_grid-1):
                grid[i,j] = .3
    return grid
#--------------------------------------------------------------------------
# Function for periodic boundaries (NOT USED)
#--------------------------------------------------------------------------
def periodic_boundary(grid):
    #
    for i in range(nghost):
        # horizontal boundary
        grid[0+i, :] = grid[l2-nghost+i+1, :]
        grid[lpy2+i, :] = grid[l1+i, :]
        # vertical boundary
        grid[:, 0+i] = grid[:, c2-nghost+i+1]
        grid[:, cpy2+i] = grid[:, c1+i]
    #
    return grid
#--------------------------------------------------------------------------
# Function for absorbing boundaries (ARTICLE DEFAULT)
#--------------------------------------------------------------------------
def absorbing_boundary(grid):
    #
    for i in range(nghost):
        # horizontal boundary
        grid[0+i, :] = 0
        grid[lpy2+i, :] = 0
        # vertical boundary
        grid[:, 0+i] = 0
        grid[:, cpy2+i] = 0
    #
    return grid
#--------------------------------------------------------------------------
# Function to update states of CA
#--------------------------------------------------------------------------
def state_update(N, I, S, R, e, c, m, v, omega, neighbors):
    # Setting mu grid
    mu = c * m * v
    # Creating grids for loop update (original ones must remain unchanged
    # and will be updated before return)
    I_up = np.array(I)
    S_up = np.array(S)
    R_up = np.array(R)
    # Updating grids
    for i in np.arange(l1, lpy2):
        for j in np.arange(c1, cpy2):
            #-------------------------------------
            # Updating I
            # Summation of Eq. 8
            I_sum = 0
            if neighbors == 4: # VonNeumann:
                for i_n, j_n in [[i-1,j], [i,j+1], [i+1,j], [i,j-1]]:
                    I_sum = I_sum + (N[i_n, j_n]/N[i,j]) * \
                            mu[i_n,j_n] * I[i_n, j_n]
            if neighbors == 8: # Moore:
                for i_n, j_n in [[i-1,j], [i-1,j+1], [i,j+1], [i+1,j+1], \
                                 [i+1,j], [i+1,j-1], [i, j-1], [i-1,j-1]]:
                    I_sum = I_sum + (N[i_n, j_n]/N[i,j]) * \
                            mu[i_n,j_n] * I[i_n, j_n]
            # Solving Eq. 8
            I_up[i,j] = (1 - e) * I[i,j] + \
                        v * S[i,j] * I[i,j] + \
                        S[i,j] * I_sum
            #-------------------------------------
            # Updating S
            # Summation of Eq. 9
            S_sum = 0
            if neighbors == 4: # VonNeumann:
                for i_n, j_n in [[i-1,j], [i,j+1], [i+1,j], [i,j-1]]:
                    S_sum = S_sum + (N[i_n, j_n]/N[i,j]) * \
                            mu[i_n,j_n] * I[i_n, j_n]
            if neighbors == 8: # Moore:
                for i_n, j_n in [[i-1,j], [i-1,j+1], [i,j+1], [i+1,j+1], \
                                 [i+1,j], [i+1,j-1], [i, j-1], [i-1,j-1]]:
                     S_sum = S_sum + (N[i_n, j_n]/N[i,j]) * \
                             mu[i_n,j_n] * I[i_n, j_n]
            # Solving Eq. 9
            S_up[i,j] = S[i,j] - \
                        omega * S[i,j] - \
                        v * S[i,j] * I[i,j] - \
                        S[i,j] * S_sum
            #-------------------------------------
            # Updating R
            # Solving Eq. 10
            R_up[i,j] = R[i,j] + \
                         e * I[i,j] + \
                         omega * S[i,j]
    #
    # Correcting I,S,R values to respect Q domain as in Eq. 6
    I_up[np.where(I_up<0)]=0.
    I_up[np.where(I_up>1)]=1.
    S_up[np.where(S_up<0)]=0.
    S_up[np.where(S_up>1)]=1.
    R_up[np.where(R_up<0)]=0.
    R_up[np.where(R_up>1)]=1.
    # Updating real grids with discretization D as in Eq. 7
    I = np.round(100*np.array(I_up))/100.
    S = np.round(100*np.array(S_up))/100.
    R = np.round(100*np.array(R_up))/100.
    #
    return I, S, R
#==========================================================================
# END OF AUXILIAR METHODS

#==========================================================================
#                           SETTING UP GRIDS
#--------------------------------------------------------------------------
# Lists for analysis
I_list = []
S_list = []
R_list = []
Omega_list=  []
#
# Creating grids
N = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
I = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
S = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
R = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
c = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
m = np.zeros([resolution + 2*nghost, resolution + 2*nghost])
#
# Applying initial Conditions
t = 0
# N
N = set_N(N, case_N)
# I
I = set_I(I, case_I)
# S
S = set_S(S, case_S)
# R
R = set_R(R, case_R) 
# e
e = e
# c
c = set_c(c, case_c)
# m
m = set_m(m, case_m)
# v
v = v
# Basic reproductive number as Eq. 16
Omega = v*np.max(S)/e
#
# Applying boundary conditions
I = absorbing_boundary(I)
S = absorbing_boundary(S)
R = absorbing_boundary(R)
N = absorbing_boundary(N)
#
# Updating lists
I_list.append(np.mean(I[l1:lpy2, c1:cpy2]))
S_list.append(np.mean(S[l1:lpy2, c1:cpy2]))
R_list.append(np.mean(R[l1:lpy2, c1:cpy2]))
Omega_list.append(Omega)
#
# First plots before entering loop
save_subplots_fig(I, \
                  I_list, S_list, R_list, \
                  t, \
                  l1, lpy2, c1, cpy2, \
                  output_folder, \
                  'grid_sir_')
save_grid_fig(I, \
                  t, \
                  l1, lpy2, c1, cpy2, \
                  output_folder, \
                  'grid')
#
# Advancing time
t = t + 1
#==========================================================================
# END OF SETTING UP GRIDS

#==========================================================================
#                   SIMULATING AN EPIDEMIC SPREADING
#--------------------------------------------------------------------------
# Auxiliar variable to get info about end of infection
# by counting for how long Omega (reprod. n.) is lower then one
Omega_counter = 0
#
# Initiating loop
while Omega_counter < 10 or t < t_max:
    #
    # Updating states
    I, S, R = state_update(N, I, S, R, e, c, m, v, omega, neighbors)
    #
    # Applying boundary conditions
    I = absorbing_boundary(I)
    S = absorbing_boundary(S)
    R = absorbing_boundary(R)
    N = absorbing_boundary(N)
    #
    # Updating reproduction number
    Omega = v*np.max(S)/e
    #
    # Updating auxiliar variable
    if Omega < 1:
        Omega_counter = Omega_counter + 1
    #
    # When to employ vaccination
    #if t==16:
    #    omega=.4
    #
    # Printing evolution 
    print('time:', t, '/ Basic repr. n.:', Omega)
    #
    # Appending results to lists 
    I_list.append(np.mean(I[l1:lpy2, c1:cpy2]))
    S_list.append(np.mean(S[l1:lpy2, c1:cpy2]))
    R_list.append(np.mean(R[l1:lpy2, c1:cpy2]))
    Omega_list.append(Omega)
    #
    # Saving outputs
    #np.savetxt(output_folder + '/I'+str(t)+'.txt', I)
    save_subplots_fig(I, \
                      I_list, S_list, R_list, \
                      t, \
                      l1, lpy2, c1, cpy2, \
                      output_folder, \
                      'grid_sir_')
    save_grid_fig(I, \
                  t, \
                  l1, lpy2, c1, cpy2, \
                  output_folder, \
                  'grid')
    #
    # Advancing time
    t = t + 1
#==========================================================================
# END OF SIMULATION

#==========================================================================
#                            SAVING OUTPUTS
#--------------------------------------------------------------------------
# Saving lists
np.savetxt(output_folder + '/I.txt', I_list)
np.savetxt(output_folder + '/S.txt', S_list)
np.savetxt(output_folder + '/R.txt', R_list)
np.savetxt(output_folder + '/Omega.txt', Omega_list)
# Final plot
sir_plot(I_list, S_list, R_list, \
           t, \
           l1, lpy2, c1, cpy2, \
           output_folder, \
           'SIR')
# Text file with input parameters
params_save(resolution, neighbors, \
          case_N, \
          case_I, case_S, case_R, \
          e, \
          case_c, \
          case_m, \
          v, \
          omega, \
          output_folder, \
          'parameters')
#==========================================================================
# END OF PROGRAM
