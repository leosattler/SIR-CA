###########################################################################
# Auxiliar file to handle outputs.
#
# INPE, Sao Jose dos Campos, SP, Brasil - July 6th, 2020
# Leonardo Sattler Cassara - leocassara@igeo.ufrj.br
###########################################################################

#==========================================================================
#                                 IMPORTS
#--------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib
import os
#==========================================================================

#==========================================================================
#                              SAVING FIGURES
#--------------------------------------------------------------------------
# Function to save snapshots of simulation
#--------------------------------------------------------------------------
def save_grid_fig(grid, \
                  t, \
                  l1, lpy2, c1, cpy2, \
                  output_folder, \
                  fig_name):
    #
    if t==0:
        os.mkdir(output_folder + '/grid_evolution')
    #
    my_cmap = matplotlib.cm.get_cmap('Greys')
    #
    plt.imshow(grid[l1:lpy2, c1:cpy2], \
               cmap=my_cmap, \
               vmin=0, vmax=1)
    #
    plt.xticks([])
    plt.yticks([])
    #
    plt.title('Day: ' + str(t))
    #
    time_str = str(t)
    while len(time_str)<4:
        time_str = '0'+time_str
    #
    plt.savefig(output_folder + '/grid_evolution/' + \
                fig_name + time_str + '.jpg', dpi=300)
    plt.close('all')
#--------------------------------------------------------------------------
# Function to save snapshots + plot of S, I and R curves
#--------------------------------------------------------------------------
def save_subplots_fig(grid, \
                      I_list, S_list, R_list, \
                      t, \
                      l1, lpy2, c1, cpy2, \
                      output_folder, \
                      fig_name):
    #
    if t==0:
        os.mkdir(output_folder + '/grid_SIR_evolution')
    #
    plt.close('all')
    #
    time_str = str(t)
    while len(time_str)<4:
        time_str = '0'+time_str
    #
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[14, 7], sharey='all')
    #
    fig.suptitle('Day: ' + str(t), size=16)
    #
    fig.subplots_adjust(wspace=0.1)
    fig.subplots_adjust(hspace=0.1)
    #
    my_cmap = matplotlib.cm.get_cmap('Greys')
    #
    ax1.imshow(grid[l1:lpy2, c1:cpy2], \
               cmap=my_cmap, \
               vmin=0, vmax=1)
    #
    ax1.set_xticks([])
    ax1.set_yticks([])
    #
    plt.xlim(0,50)
    plt.xlabel('Days since first case', size=14)
    #
    ax2_right = ax2.twinx()
    #
    ax2_right.plot(S_list, 'c-.', label='Susceptible')
    ax2_right.plot(I_list, 'r--', label='Infected')
    ax2_right.plot(R_list, 'g-', label='Recovered')
    ax2.tick_params(axis='y', which='both', \
                    labelleft=False, labelright=False)
    ax2_right.tick_params(axis='y', which='both', \
                          labelleft=False, labelright=True)
    #
    plt.ylabel('Fraction of population', size=14)
    #
    plt.ylim(0,1)
    #
    plt.savefig(output_folder + '/grid_SIR_evolution/' + \
                fig_name + time_str + '.jpg', dpi=300, bbox_inches='tight')
    plt.close('all')
#--------------------------------------------------------------------------
# Function to plot of save S, I and R curves
#--------------------------------------------------------------------------
def sir_plot(I_list, S_list, R_list, \
               t, \
               l1, lpy2, c1, cpy2, \
               output_folder, \
               fig_name):
    plt.close('all')
    #
    plt.plot(S_list, 'c-.', label='Susceptible')
    plt.plot(I_list, 'r--', label='Infected')
    plt.plot(R_list, 'g-', label='Recovered')
    #
    plt.xlim(0,50)
    plt.ylabel('Fraction of population', size=14)
    plt.xlabel('Days since first case', size=14)
    #
    plt.ylim(0,1.1)
    #
    plt.legend(loc=0)
    #
    time_str = str(t)
    while len(time_str)<4:
        time_str = '0'+time_str
    #
    plt.savefig(output_folder + '/' + \
                fig_name + '.jpg', dpi=300, bbox_inches='tight')
    plt.close('all')
#==========================================================================

#==========================================================================
#                      SAVING .txt LIST OF PARAMETERS
#--------------------------------------------------------------------------
def params_save(resolution, neighbors, \
              case_N, \
              case_I, case_S, case_R, \
              e, \
              case_c, \
              case_m, \
              v, \
              omega, \
              output_folder, \
              file_name):
    #
    file_id = output_folder + '/' + file_name + '.txt'
    #
    f = open(file_id, 'a')
    #
    f.write('GRID PARAMETERS' + '\n')
    f.write('Resolution: ' + str(resolution) + '\n')
    f.write('Number of neighbors: ' + str(neighbors) + '\n' + '\n')
    #
    f.write('INFECTION PARAMETERS' + '\n')
    # N
    if case_N == 1:
        f.write('N: ' + '100/cell' + '\n')
    if case_N == 2:
        f.write('N: ' + 'exp(j)'  + '\n')
    # I
    if case_I == 1:
        f.write('S0: ' + '.7 on center of grid'  + '\n')
    # S
    if case_S == 1:
        f.write('I0: ' + '.3 on center of grid'  + '\n')
    # R
    if case_R == 1:
        f.write('R0: ' + '0 on center of grid'  + '\n')
    # e
    f.write('e: ' + str(e) + '\n')
    # c
    if case_c == 1:
        f.write('c: ' + '1 everywhere' + '\n')
    if case_c == 2:
        f.write('c: ' + 'as article pg. 200' + '\n')
    # m
    if case_m == 1:
        f.write('m: ' + '.5 everywhere' + '\n')
    # v
    f.write('v: ' + str(v) + '\n')
    # omega
    f.write('omega: ' + str(omega) + '\n')
    #
    f.close()
#==========================================================================
