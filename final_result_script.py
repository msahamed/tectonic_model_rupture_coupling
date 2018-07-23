# example using output class in python

# required arguments are problem name and output unit name
# data directory is optional, if no argument provided assumes it is the current working directory

from __future__ import print_function
import sys, os

# dir_path = '../data/'
# sys.path.append(dir_path)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import fdfault
import pickle

# mesh parameters
maxy = 70.0
refine = 1
nby0 = 600*refine+1
nby1 = 200*refine+1

# latex parameter
# font = {
#     'family': 'serif', 
#     'serif': ['Computer Modern Roman'],
#     'weight' : 'regular',
#     'size'   : 14
#     }

# plt.rc('font', **font)
# plt.rc('text', usetex=True)
plt.style.use(['default'])
params = {
    'text.latex.preamble': ['\\usepackage{gensymb}'],
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 150,  # to adjust notebook inline plot size
    'axes.labelsize': 12, # fontsize for x and y labels (was 10)
    'axes.titlesize': 12,
    'font.size': 14, # was 10
    'legend.fontsize': 12, # was 10
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'text.usetex': False,
    'figure.facecolor':'white',
    'font.family': 'serif',
}
plt.rcParams.update(params)

color_map = 'inferno'

def initial_slip_rate(models, models_litho):
    
    # On fault slip rate + normal and shear stress
    LVfault = fdfault.output(models[0],'vfault')
    LCVfault = fdfault.output(models[1],'vfault')

    LVfault_litho = fdfault.output(models_litho[0],'vfault')
    LCVfault_litho = fdfault.output(models_litho[1],'vfault')

    LVfault.load()    
    LCVfault.load()

    LVfault_litho.load()
    LCVfault_litho.load()

    # plt.figure(figsize=(6, 4))
    plt.subplot(141)
    plt.plot(LVfault.V[0,:], LVfault.y-maxy, 'b-', label="RSN")
    plt.xlabel("Slip rate(m/s)")
    plt.ylabel("Depth (km)")
    plt.title('(a) RSN', loc='left')
    plt.xticks([0, 0.20, 0.40, 0.60])

    plt.subplot(142)
    plt.plot(LVfault_litho.V[0,:], LVfault.y-maxy, 'b-', label="RSNL")
    plt.xlabel("Slip rate(m/s)")
    plt.title('(c) RSNL', loc='left')
    plt.xticks([0, 0.20, 0.40, 0.60])
    # plt.yticks([])

    plt.subplot(143)
    plt.plot(LCVfault.V[0,:], LVfault.y-maxy, 'b-', label="RSC")
    plt.xlabel("Slip rate(m/s)")
    plt.title('(b) RSC', loc='left')
    plt.xticks([0, 0.25, 0.50, 0.75])
    # plt.yticks([])

    plt.subplot(144)
    plt.plot(LCVfault_litho.V[0,:], LVfault.y-maxy, 'b-', label="RSCL")
    plt.xlabel("Slip rate(m/s)")
    plt.title('(d) RSCL', loc='left')
    plt.xticks([0, 0.25, 0.50, 0.75])
    # plt.yticks([])
    
    # plt.tight_layout()
    plt.show()

def initial_stress_conditions(models):
    LSfault = fdfault.output(models[0],'sfault')
    LSnfault = fdfault.output(models[0],'snfault')
    SSnfault = fdfault.output(models[1],'snfault')
    Sfault = fdfault.output(models[1],'sfault')

    LSfault.load()
    Sfault.load()
    LSnfault.load()
    SSnfault.load()

    plt.figure(figsize=(5, 4))
    plt.subplot(121)
    plt.plot(LSnfault.Sn[0,:], LSnfault.y-maxy, 'b-', label="TM1")
    # plt.plot(SSnfault.Sn[0,:], SSnfault.y-maxy, 'b-', label ="TM2")
    plt.xlabel("Normal stress (MPa)")
    plt.ylabel("Depth (km)")
    # plt.legend(loc='upper left', frameon=False)
    plt.title('(a)', loc='left')
    plt.xticks([0, -100, -200, -300])

    plt.subplot(122)
    plt.plot(LSfault.S[0,:], LSfault.y-maxy, 'b-')
    # plt.plot(Sfault.S[0,:], Sfault.y-maxy, 'b-')
    plt.xlabel("Shear stress (MPa)")
    plt.title('(b)', loc='left')
    plt.xticks([0, 25, 50, 75, 100])
    plt.tight_layout()
    plt.show()

def calculate_pgv(models, models_litho):
    Lbody_vx_shallow = fdfault.output(models[0],'vxbody')
    Lbody_vy_shallow = fdfault.output(models[0],'vybody')
    Lbody_vx_deep = fdfault.output(models[1],'vxbody')
    Lbody_vy_deep = fdfault.output(models[1],'vybody')


    litho_vx_shallow = fdfault.output(models_litho[0],'vxbody')
    litho_vy_shallow = fdfault.output(models_litho[0],'vybody')
    litho_vx_deep = fdfault.output(models_litho[1],'vxbody')
    litho_vy_deep = fdfault.output(models_litho[1],'vybody')

    Lbody_vx_shallow.load()
    Lbody_vy_shallow.load()
    Lbody_vx_deep.load()
    Lbody_vy_deep.load()
    litho_vx_shallow.load()
    litho_vy_shallow.load()
    litho_vx_deep.load()
    litho_vy_deep.load()


    shallow = np.sqrt(Lbody_vx_shallow.vx[:,:, 401]**2 + Lbody_vy_shallow.vy[:,:, 401]**2)
    LV_shallow = shallow.max(axis=0)
    deep = np.sqrt(Lbody_vx_deep.vx[:,:, 401]**2 + Lbody_vy_deep.vy[:,:, 401]**2)
    LV_deep = deep.max(axis=0)

    litho_shallow = np.sqrt(litho_vx_shallow.vx[:,:, 401]**2 + litho_vy_shallow.vy[:,:, 401]**2)
    LV_litho_shallow = litho_shallow.max(axis=0)
    litho_deep = np.sqrt(litho_vx_deep.vx[:,:, 401]**2 + litho_vy_deep.vy[:,:, 401]**2)
    LV_litho_deep = litho_deep.max(axis=0)

    # load surface data
    filename = open('surface.pkl', 'rb')
    surface = pickle.load(filename)
    filename.close()

    plt.subplot(511)
    plt.plot(surface['x1'], surface['y1'], 'b-')
    plt.plot(surface['x2'], surface['y2'], 'b-')
    plt.plot(surface['xfault'], surface['yfault']-70, color='red', linestyle='--', linewidth=3)
    plt.xticks([])
    plt.ylabel('Depth (km)')
    plt.yticks([0, -5, -10])
    plt.text(55, -5, 'Footwall')
    plt.text(5, -5, 'Hanging wall')
    plt.text(0, -3, '(a)')
    plt.annotate('Fault', xy=(35, -6.5), xytext=(45, -8),
            arrowprops=dict(facecolor='black', shrink=0.005), color='red')

    plt.subplot(512)
    plt.plot(Lbody_vx_shallow.x[:, 401], LV_shallow, 'b-')
    plt.yticks([0, 1, 2])
    plt.xticks([])
    plt.text(0, 1.5, '(b) RSN')
    # plt.legend(loc='upper right', frameon=False)

    plt.subplot(513)
    plt.plot(Lbody_vx_shallow.x[:, 401], LV_litho_shallow, 'r-')
    plt.yticks([0, 10, 20])
    plt.xticks([])
    plt.text(0, 15, '(c) RSNL')
    # plt.legend(loc='upper right', frameon=False)

    plt.subplot(514)
    plt.plot(Lbody_vx_shallow.x[:, 401], LV_deep, 'b-')
    plt.yticks([0, 1, 2])
    plt.text(0, 1.1, '(d) RSC')
    plt.xticks([])
    # plt.legend(loc='upper right', frameon=False)

    plt.subplot(515)
    plt.plot(Lbody_vx_shallow.x[:, 401], LV_litho_deep, 'r-')
    plt.yticks([0, 5, 10])
    plt.xlabel('Distance across fault (km)')
    plt.text(0, 7, '(e) RSCL')
    # plt.legend(loc='upper right', frameon=False)
    
    plt.text(-14, 40, 'Peak ground velocity (m/s)', rotation=90)
    plt.show()


# ------- peak slip along the fault --------
def peak_slip(models, models_litho):
    Lfault_shallow = fdfault.output(models[0],'ufault')
    Lfault_deep = fdfault.output(models[1],'ufault')
    
    litho_shallow = fdfault.output(models_litho[0],'ufault')
    litho_deep = fdfault.output(models_litho[1],'ufault')

    normal_stress = fdfault.output(models[0],'snfault')
    shear_stress = fdfault.output(models[0],'sfault')

    Lfault_shallow.load()
    Lfault_deep.load()
    litho_shallow.load()
    litho_deep.load()
    normal_stress.load()
    shear_stress.load()

    Lfault_shallow_slip = Lfault_shallow.U.max(axis=0)
    Lfault_deep_slip = Lfault_deep.U.max(axis=0)
    
    litho_shallow_slip = litho_shallow.U.max(axis=0)
    litho_deep_slip = litho_deep.U.max(axis=0)
    depths = Lfault_shallow.y

    # plt.figure(figsize=(6.20, 4))
    plt.subplot(131)
    plt.plot(Lfault_shallow_slip, depths-maxy, 'b-', label="RSN")
    plt.plot(litho_shallow_slip, depths-maxy, 'r-', label="RSNL")
    plt.xlabel('Peak Slip (m)')
    plt.ylabel('Depth (km)')
    plt.legend(loc='lower right', frameon=False)
    plt.title('(a)', loc='left')

    plt.subplot(132)
    plt.plot(Lfault_deep_slip, depths-maxy, 'b-', label="RSC")
    plt.plot(litho_deep_slip, depths-maxy, 'r-', label="RSCL")
    plt.xlabel('Peak Slip (m)')
    plt.legend(loc='lower right', frameon=False)
    plt.title('(b)', loc='left')

    plt.subplot(133)
    plt.plot(np.absolute(shear_stress.S[0,:])/np.absolute(normal_stress.Sn[0,:]), depths-maxy, 'b-', label="T1C")
    plt.xlabel('Shear stress/Normal stress')
    # plt.xticks([0, 100, 200])
    plt.title('(c)', loc='left')
    plt.tight_layout()
    plt.show()

# ------- spatio_temporal_slip_rate -------
def spatio_temporal_slip_rate(models, models_litho):

    Lfault_shallow = fdfault.output(models[0],'vfault')
    Lfault_deep = fdfault.output(models[1],'vfault')

    litho_shallow = fdfault.output(models_litho[0],'vfault')
    litho_deep = fdfault.output(models_litho[1],'vfault')

    Lfault_shallow.load()
    Lfault_deep.load()
    litho_shallow.load()
    litho_deep.load()

    times = Lfault_shallow.t
    depths = Lfault_shallow.y
    time_len = int(len(times))
    meshX , meshY = np.meshgrid(times[:time_len], depths)

    fig = plt.figure(figsize=(5.5, 4))
    plt.subplot(141)
    plt.pcolor(meshX , meshY-maxy, (Lfault_shallow.V).T, vmax = 10, vmin = 0, cmap = color_map)
    plt.xlabel('Time (s)')
    plt.xticks([0,1, 2, 3, 4])
    plt.title('(a) RSN', loc='left')
    plt.ylabel('Depth (km)')
    # cb = plt.colorbar()

    plt.subplot(142)
    plt.pcolor(meshX , meshY-maxy, (litho_shallow.V).T, vmax = 10, vmin = 0, cmap = color_map)
    plt.xlabel('Time (s)')
    plt.yticks([])
    plt.xticks([0,1, 2, 3, 4])
    plt.title('(b) RSNL', loc='left')

    plt.subplot(143)
    ax1 = plt.pcolor(meshX , meshY-maxy, (Lfault_deep.V).T, vmax = 10, vmin = 0,cmap = color_map)
    plt.xlabel('Time (s)')
    plt.yticks([])
    plt.xticks([0,1, 2, 3, 4])
    plt.title('(c) RSC', loc='left')
    # cb = plt.colorbar()

    plt.subplot(144)
    ax1 = plt.pcolor(meshX , meshY-maxy, (litho_deep.V).T, vmax = 10, vmin = 0, cmap = color_map)
    plt.xlabel('Time (s)')
    plt.yticks([])
    plt.xticks([0,1, 2, 3, 4])
    plt.title('(d) RSCL', loc='left')
    
    cb = fig.colorbar(ax1, ticks=[0, 5, 10])
    cb.set_label('Slip veolocity (m/s)')
    plt.tight_layout()
    plt.show()

# ------ sliprate in different Times -----
def plot_slip_rate(models):
    steps = [0, 2, 4, 6, 8, 10]
    c = ['b', 'g', 'r', 'c', 'm', 'k', 'y']

    LVfault = fdfault.output(models[0],'vfault')
    SVfault = fdfault.output(models[1],'vfault')
    LVfault_coh = fdfault.output(models[2],'vfault')
    SVfault_coh = fdfault.output(models[3],'vfault')
    LVfault.load()
    SVfault.load()
    LVfault_coh.load()
    SVfault_coh.load()

    plt.figure()
    plt.subplot(1,2,1)
    for i in range(len(steps)):
        c1 = str(c[i]+'-')
        c2 = str(c[i]+'--')

        plt.plot(LVfault.V[steps[i],:], LVfault.y-maxy, c1, label= str(round(LVfault.t[steps[i]],0))+ ' s')
        plt.plot(SVfault.V[steps[i],:], SVfault.y-maxy, c2)

    plt.title('(a)', loc='left')
    plt.ylabel("Distance along fault (km)")
    plt.xlabel("Slip rate (m/s)")
    plt.legend(loc="upper right", frameon=False)
    plt.grid(True)
    
    plt.subplot(1,2,2)
    for i in range(len(steps)):
        c1 = str(c[i]+'-')
        c2 = str(c[i]+'--')

        plt.plot(LVfault_coh.V[steps[i],:], LVfault_coh.y-maxy, c1, label= str(round(LVfault_coh.t[steps[i]],0))+ ' s')
        plt.plot(SVfault_coh.V[steps[i],:], SVfault_coh.y-maxy, c2)

    plt.title('(b)', loc='left')
    # plt.ylabel("Depth")
    plt.xlabel("Slip rate (m/s)")
    plt.legend(loc="upper right", frameon=False)
    plt.grid(True)
    plt.show()


# color particle velocity
def set_color_bar(a, lbl, fig1):
    fig1.subplots_adjust(hspace=0.5)
    axes =[0.91, 0.1, 0.02, 0.8]
    cbar_ax = fig1.add_axes(axes)
    cb = fig1.colorbar(a, orientation='vertical', cax=cbar_ax)
    tick_locator = ticker.MaxNLocator(nbins=4)
    cb.locator = tick_locator
    color_bar_label = lbl+' particle velocity($m/s$)'
    cb.set_label(color_bar_label)
    cb.update_ticks()

#--------------  Horizontal particle velocity -------------
def print_horizontal_velocity(models, steps):
    row = len(steps)
    col = 2
    mx = 0.5
    mn = mx*(-1)
    Lbody_vx = fdfault.output(models[0],'vxbody')
    Sbody_vx = fdfault.output(models[1],'vxbody')
    Lbody_vx.load()
    Sbody_vx.load()
    # maxy = max(Lbody_vx.y)
    fig1 = plt.figure(figsize=(9, 5))
    indx = 0
    label = 'Horizontal'
    
    for i in range(row):
        # long term model
        indx = indx + 1
        ax1 = fig1.add_subplot(row, col, indx)
        a = ax1.pcolor(Lbody_vx.x, Lbody_vx.y, Lbody_vx.vx[steps[i], :, :], cmap=color_map, vmax= mx, vmin=mn)
        ax1.set_yticklabels([-10, 0])
        # ax1.tick_params(axis='x', which='both', bottom='off', top='off',right='off')
        # ax1.set_title('Time: '+ str(round(Lbody_vx.t[steps[i]],0))+'  s', loc='left')
        
        if indx == 1:
            ax1.set_title('(a) \t\t\t'+ 't = '+str(round(Lbody_vx.t[steps[i]],2))+'  s', loc='left')
            ax1.set_xticklabels([])
            ax1.set_ylabel('Depth (km)')
        if indx == (row*2-1):
            ax1.set_xlabel('Distance (km)')
            ax1.set_ylabel('Depth (km)')
            ax1.set_title('t = '+str(round(Lbody_vx.t[steps[i]],2))+'  s', loc='left')
        elif indx > 1 and indx<(row*2-1):
            ax1.set_xticklabels([])
            ax1.set_title('t = '+str(round(Lbody_vx.t[steps[i]],2))+'  s', loc='left')


        # short term model
        indx = indx + 1
        ax1 = fig1.add_subplot(row, col, indx)
        a = ax1.pcolor(Sbody_vx.x, Sbody_vx.y, Sbody_vx.vx[steps[i], :, :], cmap=color_map, vmax=mx, vmin=mn)
        ax1.set_yticklabels([])
        # ax1.tick_params(axis='x', which='both', bottom='off', top='off', right='off')
        
        if indx == row*2:
            ax1.set_xlabel('Distance (km)')
        elif indx == 2:
            ax1.set_title('(b)', loc='left')
        else:
            ax1.set_xticklabels([])

    # fig1.savefig('rupture_vx_without_coh.png')   # save the figure to file
    # plt.close(fig1)
    set_color_bar(a, label, fig1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
    plt.show()

#--------------  Horizontal particle velocity -------------
def print_vertical_velocity(models, steps, color_map):
    # steps = [50, 100, 150, 200]
    row = len(steps)
    col = 2
    mx = 0.5
    mn = mx*(-1)
    Lbody_vy = fdfault.output(models[0],'vybody')
    Sbody_vy = fdfault.output(models[1],'vybody')
    Lbody_vy.load()
    Sbody_vy.load()
    fig1 = plt.figure()
    indx = 1
    label = 'Vertical'
    
    for i in range(row):
        print(indx)
        # long term model
        ax1 = fig1.add_subplot(row, col, indx)
        a = ax1.pcolor(Lbody_vy.x, Lbody_vy.y, Lbody_vy.vy[steps[i], :, :], cmap=color_map, vmax= mx, vmin=mn)
        ax1.set_yticklabels([-10, '', '', '', '', 0])
        ax1.tick_params(axis='x', # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        right='off')
        ax1.set_title('Time: '+ str(round(Lbody_vy.t[steps[i]],0))+'  s', loc='left')
        ax1.set_ylabel('Depth (km)')
        if indx==7:
            set_color_bar(a, label, fig1)
            ax1.set_xlabel('Distance (km)')
        elif indx == 1:
            ax1.set_title('(a) \n Time: '+ str(round(Lbody_vy.t[steps[i]],0))+'  s', loc='left')
        else:
            ax1.set_xticklabels([])

        # short term model
        ax1 = fig1.add_subplot(row, col, indx+1)
        a = ax1.pcolor(Sbody_vy.x, Sbody_vy.y, Sbody_vy.vy[steps[i], :, :], cmap=color_map, vmax=mx, vmin=mn)
        ax1.set_yticklabels([])
        ax1.tick_params(axis='x', # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',        # ticks along the top edge are off
        right='off')
        if indx==7:
            set_color_bar(a, label, fig1)
            ax1.set_xlabel('Distance (km)')
        elif indx ==1:
            ax1.set_title('(b) \n', loc='left')
        else:
            ax1.set_xticklabels([])
        indx+=2

    # fig1.savefig('rupture_vy_without_coh.png')   # save the figure to file
    # plt.close(fig1)
    plt.show()


def main():
    steps =  [7, 15, 20, 40]
    with_cohesion = True
    models = ['T1N_68125', 'T1C_68125']
    # model_no_coh = ['T1N_68125', 'T2N_68125']
    models_litho = ['longterm_litho_shallow', 'longterm_litho_deep']

    # plot_slip_rate(slip_models)
    initial_stress_conditions(models, models_litho)
    # initial_slip_rate(models, models_litho)
    # calculate_pgv(models, models_litho)
    # peak_slip(models, models_litho)
    # spatio_temporal_slip_rate(models, models_litho)
    # print_horizontal_velocity(models,steps)
    # print_vertical_velocity(models, steps)
    # make_horizontal_animation()

if __name__ == "__main__":
    main()
