import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
import copy
import os
import warnings
warnings.filterwarnings("ignore", message="Multiprocessing-backed parallel loops cannot be nested")
from spectrometer.Experiment import MeV2m0s, vel2gamma
from spectrometer.Experiment import Experiment
import spectrometer.plot_result as pr
from spectrometer.geometry import *
from spectrometer.physics import velocity_devider


def _log_scale_vectors(Fx, Fy, Fz, len_min=0.3, len_max=1.0, eps=1e2):
    F_mag = np.sqrt(Fx**2 + Fy**2 + Fz**2)
    F_mag_safe = np.maximum(F_mag, eps)

    log_mag = np.log10(F_mag_safe)
    log_min = np.log10(eps)
    log_max = np.max(log_mag)

    if log_max > log_min:
        scaled_len = len_min + (len_max - len_min) * (log_mag - log_min) / (log_max - log_min)
    else:
        scaled_len = np.full_like(log_mag, len_min)

    Fx_s = (Fx / F_mag_safe) * scaled_len
    Fy_s = (Fy / F_mag_safe) * scaled_len
    Fz_s = (Fz / F_mag_safe) * scaled_len
    return Fx_s, Fy_s, Fz_s



def plot_field_quiver(ax,integrator):
    Bx0_T     = integrator.Bx0_T
    Ex0_Vm    = integrator.Ex0_Vm 
    width_mm  = integrator.width_mm
    depth_mm  = integrator.depth_mm
    height_mm = integrator.height_mm
    yoke      = integrator.yoke
    fringe    = integrator.fringe
    shield_mm = integrator.shield_mm


    def _eps_from_scale(x0):
        x0 = abs(x0)
        if x0 == 0:
            return 1.0
        return 10**abs(round(np.log10(x0)))


    field_vec = []
    if Ex0_Vm !=0:
        rutin = 2
    else:
        rutin = 1
    field = integrator
        
        
    quive_len_vec = [4,4]
    shift_vec = [0,2]
    
    color_vec = ["magenta", "cyan", ]
    for r in range(rutin):
        
        quive_len = quive_len_vec[r]
        color = color_vec[r]
        shift = shift_vec[r]
    # Gread resolution
        nx, ny, nz = 20, 50, 5
        
    # Gread generation
        X_width = np.linspace(-width_mm/2+shift, width_mm/2 - quive_len+shift, nx) * 1e-3
        if yoke:
            if fringe:
                Y_depth = np.linspace(min(-max(shield_mm,15),field.R0_mm[1]), depth_mm * 1, ny) * 1e-3 
            else:
                Y_depth = np.linspace(0, depth_mm * 1, ny) * 1e-3 
        else:
            if fringe:
                Y_depth = np.linspace(min (-max(shield_mm,15),field.R0_mm[1]), depth_mm + max(shield_mm,15,abs(field.R0_mm[1])), ny) * 1e-3 
            else:
                Y_depth = np.linspace(0, depth_mm, ny) * 1e-3 
                
        Z_height = np.linspace(-height_mm/2, height_mm/2, nz) * 1e-3
        
        X, Y, Z = np.meshgrid(X_width, Y_depth, Z_height, indexing='ij')
        
        Fx = np.zeros_like(X)
        Fy = np.zeros_like(X)
        Fz = np.zeros_like(X)
        
    
    # Caclulating the magnetic field in every point
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    R_current_m = np.array([X[i,j,k], Y[i,j,k], Z[i,j,k]])
                    if r == 0:
                        F_vec = field._get_magnetic_field(R_current_m)
                    elif r==1:
                        F_vec = field._get_electric_field(R_current_m)
                    Fx[i,j,k] = F_vec[0]
                    Fy[i,j,k] = F_vec[1]
                    Fz[i,j,k] = F_vec[2]
                    
        X_mm = X * 1e3
        Y_mm = Y * 1e3
        Z_mm = Z * 1e3
        
    # plot magentic filed arrows
    
        # norm = np.sqrt(Fx**2 + Fy**2 +Fz**2) 
        # Fx = Fx/norm
        # Fy = Fy/norm
        # Fz = Fz/norm

        if r == 0:
            eps = _eps_from_scale(integrator.Bx0_T)
        else:
            eps = _eps_from_scale(integrator.Ex0_Vm)
        Fx, Fy, Fz = _log_scale_vectors(Fx, Fy, Fz, eps=eps)
        ax.quiver(X_mm, Y_mm, Z_mm, 
                  Fx, Fy, Fz, 
                  length=quive_len,       
                  normalize = False,
                  color=color, 
                  alpha=0.6, 
                  arrow_length_ratio = 1)
    
def one_trajectory(integrator,integator_name, experiment, engs_MeV):
    loc_integrator = copy.deepcopy(integrator)
    loc_experiment = copy.deepcopy(experiment)
    
    loc_experiment.q_eng_MeV = engs_MeV
    loc_experiment._evaluate_exp_paramas()
    params = loc_experiment.params
    loc_integrator._update_params(params)
    
    integrator._update_params(params)
    if integator_name == "Euler" or integator_name == ["Euler"]:
        R_mm,v_m0s,gamma= loc_integrator._euler()
    elif integator_name == "RK2" or integator_name == ["RK2"]:
            R_mm,v_m0s,gamma= loc_integrator._RK2()
    elif integator_name == "RK4_Lin" or integator_name == ["RK4_Lin"]:
            R_mm,v_m0s,gamma= loc_integrator._RK4_Linear()
    elif integator_name == "RK4_Herm" or integator_name == ["RK4_Herm"]:
            R_mm,v_m0s,gamma= loc_integrator._RK4_Hermit()      
    elif integator_name == "Boris" or integator_name == ["Boris"]:
            R_mm,v_m0s,gamma= loc_integrator._Boris_pusher()
    elif integator_name == "Boris_Coll" or integator_name == ["Boris_Coll"]:
            R_mm,v_m0s,gamma= loc_integrator._Collocated_Boris_pusher()
            
    return R_mm,v_m0s,gamma

def basic_simulator(integrator,integator_name, experiment, engs_MeV,show_tragectories):
        
        
        
    if show_tragectories == True:
        
        spec = Spectrometer_Budy(experiment.height_mm, experiment.width_mm, experiment.depth_mm,
                                 experiment.pinhole_dia_mm, experiment.shield_mm, experiment.yoke, experiment.shield_mm)
        fig,ax = spec._draw_spec()
        if  True:
            plot_field_quiver(ax,integrator)
    experiment.solution = integrator.solution
    
    if len(engs_MeV.shape)==1:
        rutin = 1
        
    if len(engs_MeV.shape)==2:
        rutin = engs_MeV.shape[0]
        
    n_cpus = os.cpu_count();
    n_jobs = n_cpus -1 if rutin > 1 else 1
    
    engs_mat_MeV = np.array([eng for eng in engs_MeV ])
    engs_mat_MeV = np.atleast_2d(engs_mat_MeV)
    R_vec_mm = []; v_vec_m0s = []; gamma_vec = []; label = []
    

    results = Parallel(n_jobs=n_jobs, backend="multiprocessing")(
        delayed(one_trajectory)(integrator,integator_name, experiment, engs_MeV) for engs_MeV in engs_mat_MeV )
    
    R_vec_mm  = [r[0] for r in results]
    v_vec_m0s = [r[1] for r in results]
    gamma_vec = [r[2] for r in results]

        
    if show_tragectories == True:
        pos = ax.get_position()

        for r in range(len(R_vec_mm)):
            R_mm = R_vec_mm[r]
            label = str(np.round(np.linalg.norm(engs_mat_MeV[r]),2))
            pr.trajectory_plot(ax, R_mm,label)
            print(R_mm[-1])
    
        ax.legend(title='Particle total energy',title_fontsize=25, loc='center left',
                  bbox_to_anchor=(1.2 , 0.5), fontsize=18)

        
        R_vec_mm, v_vec_m0s, gamma_vec, ax
    return R_vec_mm, v_vec_m0s, gamma_vec      
        
def energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, engs_MeV,N_steps_vec):
    
    params = experiment.params
    fig = plt.figure()
    
    manager = plt.get_current_fig_manager()
    try:
        window = manager.window
     
        window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
    except Exception as e:
        print("Could not resize window:", e)

    ax = fig.add_subplot(111)
    
    colors = {"Boris": "magenta", "Boris_Coll": "red","RK2":"blue", "RK4_Lin": "green", "RK4_Herm": "Olive", "Euler": "cyan"}
    for i in range (len(integrators_name)):
        integrator_name = integrators_name[i]
        R_vec_mm = [];  v_vec_m0s = []; 
        gamma_vec = []; gamma_retio_vec = []; gamma_end_vec = [];
        
        for n in range(len(N_steps_vec)):
            # print (N_steps_vec[n])
            params["N_steps"] = N_steps_vec[n]
            params = experiment._update_params(params)
            integrator._update_params(params)
            R_engs_mm, v_engs_m0s, gamma_engs_vec = basic_simulator (integrator,integrator_name, experiment, engs_MeV,show_tragectories=False)
            gamma_retio = []
            norm_engs_MeV = []
            gamma0 = []
            gamma_end = []
            for g in range(len(gamma_engs_vec)):
                # print(gamma_engs_vec[n])
                norm_engs_MeV.append(np.linalg.norm(engs_MeV[g]))
                gamma_retio.append(gamma_engs_vec[g][-1]/gamma_engs_vec[g][0])
                gamma_end.append(gamma_engs_vec[g][-1])
                gamma0.append(gamma_engs_vec[g][0])
            
            R_vec_mm.append(R_engs_mm);
            v_vec_m0s.append(v_engs_m0s);
            gamma_vec.append(gamma_engs_vec)
            gamma_retio_vec.append(gamma_retio)
            gamma_end_vec.append(gamma_end)
            norm_engs_MeV = np.array(norm_engs_MeV)
            
        gamma_retio_vec = np.array(gamma_retio_vec)
        gamma_end_vec = np.array(gamma_end_vec)
         
        # color = colors[i]
        # ax.plot(N_steps_vec, gamma_MeV,color=color)    
        color = colors[integrator_name]
        for g in range(gamma_end_vec.shape[1]):
            
            label = integrator_name if g == 0 else None   # legend entry one for each integrator
            ax.scatter(N_steps_vec, gamma_end_vec[:, g], color=color, alpha=0.6, s=300,
                       edgecolor='black', linewidth=0.5, label=label)
            if i==0:
                for n in range(len(N_steps_vec)):
                    label_0="$\gamma(0)$"  if n==0 and g==0 else None
                    ax.scatter(N_steps_vec[n], gamma0[g], color="black", marker = "+", alpha=1, s=1000, 
                               linewidth=0.5, label=label_0)
                
                ##### NEED TO CHECK - Somthing is weared
    
    
    ax.set_xlabel(r'$N_{steps}$',fontsize=25)
    ax.set_ylabel(r'$\gamma(exit)$',fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=18)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    
    ax.legend(title='Integrator',title_fontsize=25,fontsize=18)

def convergence_test (integrator,integrators_name, experiment, engs_MeV,N_steps_vec):
    params = experiment.params
    
    
    
    colors = {"Boris": "magenta", "Boris_Coll": "red","RK2":"blue", "RK4_Lin": "green", "RK4_Herm": "Olive", "Euler": "cyan"}
    R_ends_mm = []
    for i in range (len(integrators_name)):
        integrator_name = integrators_name[i]
        R_vec_mm = [];  v_vec_m0s = []; 
        gamma_vec = []; gamma_retio_vec = []; gamma_end_vec = [];
        
        for n in range(len(N_steps_vec)):
            # print (N_steps_vec[n])
            params["N_steps"] = N_steps_vec[n]
            params = experiment._update_params(params)
            integrator._update_params(params)
            R_engs_mm, v_engs_m0s, gamma_engs_vec = basic_simulator (integrator,integrator_name, experiment, engs_MeV,show_tragectories=False)
    
            R_vec_mm.append(R_engs_mm)
            v_vec_m0s.append(v_engs_m0s)
            gamma_vec.append(gamma_engs_vec)
        
        R_end_mm =  [r[-1][-1] for r in R_vec_mm]
        R_ends_mm.append(R_end_mm)
    
    R_ends_mm = np.array(R_ends_mm)
    
    
    

    
    
    fig = plt.figure()    
    manager = plt.get_current_fig_manager()
    try:
        window = manager.window
     
        window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
    except Exception as e:
        print("Could not resize window:", e)

    ax = fig.add_subplot(111)
    def _power_law(x, a, b):
        return a * np.power(x, b)
    
    def _linera_rig(x,a,b):
        return a * x + b
    Y_analitic = integrator._analitic_sol_vel2dist()
    n_fine = np.linspace(N_steps_vec[0],N_steps_vec[-1],1001)
    power_map = {"Boris": -1,"Boris_Coll": -2 ,"RK2": -2, "RK4_Lin" : -2 ,"RK4_Herm": -4, "Euler": -1}
    count = 0
    for n in range (len(integrators_name)):
        count += 0.05
        integrator_name = integrators_name[n]
        # Y_min = R_ends_mm[n][-1][1]
        # epsilon = np.abs([Y_min - r[1] for r in R_ends_mm[n]])
        epsilon = np.abs([Y_analitic - r[1] for r in R_ends_mm[n]])
        # print(f"{integrators_name[n]}: epsilon = {epsilon}")
        power = power_map[integrator_name]
        log_eps = np.log10(epsilon[:-1])
        log_N = np.log10(N_steps_vec[:-1])
        popt, pcov = curve_fit(_linera_rig, log_N, log_eps, p0=[power,0],bounds=([-6, -np.inf],[  0,  np.inf]),maxfev=10000)
        a_fit, b_fit = popt
        print(a_fit)
        # eps_fine = _power_law(n_fine,a_fit,b_fit)
        eps_fine = 10**(_linera_rig(np.log10(n_fine),a_fit,b_fit))
        
        color = colors[integrator_name]
        N_mid_log = 10**((0.5+count)*(np.log10(N_steps_vec[0]) + np.log10(N_steps_vec[-2])))
        eps_mid_on_fit = 10**(_linera_rig(np.log10(N_mid_log), a_fit, b_fit))
        
        
        ax.text(N_mid_log , eps_mid_on_fit*10,
                rf"$p = {a_fit:.2f}$",
                color=color, fontsize=25, ha='center', va='bottom',
                fontweight='bold')
        label = integrator_name
        ax.scatter(N_steps_vec[:], epsilon[:], color=color, alpha=0.6, s=300,
                   edgecolor='black', linewidth=0.5, label=label)
        ax.plot(n_fine,eps_fine, color=color, linewidth=3)
        
        ax.set_xlabel(r'$N_{steps}$',fontsize=25)
        ax.set_ylabel(r'$\epsilon$',fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.tick_params(axis='both', which='minor', labelsize=18)
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        ax.legend(title='Integrator',title_fontsize=25,fontsize=18)
    
    
    return    
def Non_collimataed_beams(integrator,integrator_name, experiment,N_steps, total_energy,number_of_experiments, angular_distribution, angular_scale = 0.01,show_tragectories = True, show_map = True):
    params = experiment.params
    Noe = number_of_experiments
    
    if angular_distribution == "Exponential":
        
        phi_rad = np.random.exponential(np.pi * angular_scale,Noe) * np.random.choice([-1, 1], size=Noe)  
        theta_rad = np.random.exponential(np.pi/2 *  angular_scale,Noe) * np.random.choice([-1, 1], size=Noe)  
        # total_energy_vec = total_energy * np.ones(Noe)
        
        
    elif angular_distribution == "Gaussian":
        phi_rad = np.random.normal(loc=0.0, scale=np.pi*angular_scale, size=Noe)
        theta_rad = np.random.normal(loc=0.0, scale=np.pi/2*angular_scale, size=Noe)
    
    elif angular_distribution == "Uniform":
        phi_rad = np.random.uniform(low=-np.pi*angular_scale, high = np.pi*angular_scale, size=Noe)
        theta_rad = np.random.uniform(low=-np.pi/2*angular_scale, high=-np.pi/2*angular_scale, size=Noe)
        
    elif angular_distribution == "None":
        phi_rad = np.zeros(Noe)
        theta_rad = np.zeros(Noe)
        
    engs_MeV = velocity_devider (total_energy,theta_rad ,phi_rad)
    R_engs_vec_mm = []
    v_engs_vec_vec_m0s = []
    gamma_engs_vec = []
    # for n in range(Noe):    
        # params['q_eng_MeV'] = engs_MeV[n]
        # params = experiment._update_params(params)
        # integrator._update_params(params)
        
    R_engs_mm, v_engs_m0s, gamma_engs = basic_simulator (integrator,integrator_name, experiment, engs_MeV,show_tragectories=show_tragectories)
     
    R_engs_vec_mm.append(R_engs_mm)
    v_engs_vec_vec_m0s.append(v_engs_m0s)
    gamma_engs_vec.append(gamma_engs)
    # pr.trajectory_plot (ax, R_vec_mm,lable = None)  
    if show_map == True:
        fig = plt.figure()
        
        manager = plt.get_current_fig_manager()
        try:
            window = manager.window
         
            window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
        except Exception as e:
            print("Could not resize window:", e)
    
        ax = fig.add_subplot(111)
        R_engs_im_mm = np.zeros([Noe,3])
        count = 0
        x_exit_vec_mm = []
        y_exit_vec_mm = []
        for n in range(Noe):
            R_engs_im_mm[n] = R_engs_mm[n][-1]
            if  abs( abs(R_engs_im_mm[n][2])-experiment.height_mm/2) < 1e-2:
                count +=1
                # x_exit_vec_mm.append(R_engs_im_mm[n,0])
                # y_exit_vec_mm.append(R_engs_im_mm[n,1])
        x_exit_vec_mm = R_engs_im_mm[:,0]
        y_exit_vec_mm = R_engs_im_mm[:,1]
        # x_exit_vec_mm = np.array(x_exit_vec_mm
        # y_exit_vec_mm = np.array(y_exit_vec_mm)
        pecent = count/Noe * 100 
            
        
        
        
        plt.hist2d(x_exit_vec_mm, y_exit_vec_mm, bins=100,
               range=[[-experiment.width_mm/2, experiment.width_mm/2], [0, experiment.depth_mm]],
               cmap='viridis')
        
        plt.colorbar(label='counts')
        plt.xlabel('Width (mm)',fontsize=25)
        plt.ylabel('Depth (mm)',fontsize=25)
        ax.tick_params(axis='both', labelsize=18)
        ax.axis('equal')
        ax.set_title(f"{pecent} of the particle passed the ditector with {angular_distribution} angular distribution havig scale of {angular_scale}*$\pi$ for $\phi$ and {angular_scale}*$\pi/2$ for $\theta$")
    

    return    R_engs_mm     

def real_beams(integrator,integrator_name, experiment,N_steps, total_energy_scalar_MeV,
               number_of_experiments,energy_distribution,energy_scale, angular_distribution, angular_scale = 0.01,show_tragectories = True, show_map = True):
    Noe = number_of_experiments
    if energy_distribution == "Exponential":
        total_energys = np.random.exponential(total_energy_scalar_MeV*energy_scale, Noe)
    
    elif energy_distribution == "Gaussian":
        total_energys = np.random.normal(loc=total_energy_scalar_MeV, scale=total_energy_scalar_MeV * energy_scale, size=Noe)
    
    elif energy_distribution == "Uniform":
       total_energys = np.random.uniform(low=0, high = total_energy_scalar_MeV*2, size=Noe)
      
    elif energy_distribution == "DLA":
         total_energys = energies_MeV = sc.stats.gamma.rvs(a=1.5, scale=total_energy_scalar_MeV, size=Noe)   
       
    elif energy_distribution == "None":
       total_energys = total_energy_scalar_MeV * np.ones(Noe)
    
    
    number_of_angular_experiments = 1
    R_engs_mat_mm =[]
    
    n_cpus = os.cpu_count();
    n_jobs = n_cpus
    total_energys = np.atleast_2d(total_energys)
    results = Parallel(n_jobs=n_jobs, backend="multiprocessing")(
        delayed(Non_collimataed_beams)(integrator,integrator_name, experiment,N_steps, total_energy,
                                       number_of_angular_experiments, angular_distribution, angular_scale = angular_scale,
                                       show_tragectories = False, show_map = False)  for total_energy in total_energys) 
    
    R_engs_mat_mm  = [r for r in results[0]]
    
    
    
    # for te in range(Noe):
    #     total_energy = total_energys[te]
    #     R_engs_mm = Non_collimataed_beams(integrator,integrator_name, experiment,N_steps, total_energy,number_of_angular_experiments, angular_distribution, angular_scale = angular_scale,show_tragectories = False, show_map = False)
    #     R_engs_mat_mm.append(R_engs_mm)
        
    
    R_engs_im_mm = np.zeros([Noe,3])
    count = 0
    x_exit_vec_mm = []
    y_exit_vec_mm = []
    for n in range(Noe):
        R_engs_im_mm[n] = R_engs_mat_mm[n][-1]
        if  abs( abs(R_engs_im_mm[n][2])-experiment.height_mm/2) < 1e-2:
            count +=1
            # x_exit_vec_mm.append(R_engs_im_mm[n,0])
            # y_exit_vec_mm.append(R_engs_im_mm[n,1])
    x_exit_vec_mm = R_engs_im_mm[:,0]
    y_exit_vec_mm = R_engs_im_mm[:,1]
    # x_exit_vec_mm = np.array(x_exit_vec_mm
    # y_exit_vec_mm = np.array(y_exit_vec_mm)
    pecent = count/Noe * 100 
            
        # if show_map==True:
    if show_map == True:
        fig = plt.figure()
        
        manager = plt.get_current_fig_manager()
        try:
            window = manager.window
         
            window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
        except Exception as e:
            print("Could not resize window:", e)
    
        ax = fig.add_subplot(111)  
        
        h = plt.hist2d(x_exit_vec_mm, y_exit_vec_mm, bins=100,
               range=[[-experiment.width_mm/2, experiment.width_mm/2], [0, experiment.depth_mm]],
               cmap='viridis')
        
        cbar = plt.colorbar(h[3], ax=ax)
        cbar.set_label('counts', fontsize=25)   
        cbar.ax.tick_params(labelsize=18)
        plt.xlabel('Width (mm)',fontsize=25)
        plt.ylabel('Depth (mm)',fontsize=25)
        ax.tick_params(axis='both', labelsize=18)
        ax.axis('equal')
        ax.set_title(f"{np.round(pecent,2)} of the particle passed the ditector with {angular_distribution} angular distribution havig scale of {angular_scale}*$\pi$ for $\phi$ and {angular_scale}*$\pi/2$ for $\theta$")
    
    # R_engs_return_mm = []
    # count = 0
    # for n in range(Noe):
    #     R_engs_return_mm.append(R_engs_mat_mm[n][0][-1])
    #     if  abs( abs(R_engs_return_mm[n][2])-experiment.height_mm/2) < 1e-2:
    #         count +=1
    #         # x_exit_vec_mm.append(R_engs_im_mm[n,0])
    #         # y_exit_vec_mm.append(R_engs_im_mm[n,1])
    # pecent = count/Noe * 100 
    return np.array(R_engs_im_mm), pecent
        
        
def TNSA(integrator,integrator_name, experiment,N_steps, total_energy_scalars_MeV,  mq_vec_kg , q_vec_C,
               number_of_experiments,energy_distribution,energy_scale, angular_distribution, angular_scale = 0.01,show_tragectories = True, show_map = True):
    params = experiment.params
    # mq_vec_kg = np.array
    
    R_engs_im_mat_mm=[]
    pecents_vec=[]
    for i in range(len(mq_vec_kg)):
        params["q_C"] = q_vec_C[i]
        params["m_kg"] = mq_vec_kg[i]
        total_energy_scalar_MeV = total_energy_scalars_MeV[i]
        params = experiment._update_params(params)
        integrator._update_params(params)
        
        R_engs_im_mm, pecent = real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
                   N_steps = N_steps, total_energy_scalar_MeV = total_energy_scalar_MeV, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = False)
        
        R_engs_im_mat_mm.append(R_engs_im_mm)
        pecents_vec.append(pecent)
        
        
    R_engs_im_mat_mm = np.array(R_engs_im_mat_mm)
    if len(R_engs_im_mat_mm.shape) == 3:
        x_exit_vec_mm = R_engs_im_mat_mm[:, :, 0].flatten()
        y_exit_vec_mm = R_engs_im_mat_mm[:, :, 1].flatten()
    else:
        x_exit_vec_mm = R_engs_im_mat_mm[:, 0]
        y_exit_vec_mm = R_engs_im_mat_mm[:, 1]
        
        
    if show_map == True:
        fig = plt.figure()
        
        manager = plt.get_current_fig_manager()
        try:
            window = manager.window
         
            window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
        except Exception as e:
            print("Could not resize window:", e)
        ax = fig.add_subplot(111)
        
        h = plt.hist2d(x_exit_vec_mm,y_exit_vec_mm, bins=100,
               range=[[-experiment.width_mm/2, experiment.width_mm/2], [0, experiment.depth_mm]],
               cmap='viridis')
        
        cbar = plt.colorbar(h[3], ax=ax)
        cbar.set_label('counts', fontsize=25)   
        cbar.ax.tick_params(labelsize=18)
        plt.ylabel('Dept (mm)',fontsize=25)
        plt.xlabel('Widthh (mm)',fontsize=25)
        ax.tick_params(axis='both', labelsize=18)
        ax.axis('equal')
        ax.set_title(f"{pecent} of the particle passed the ditector with {angular_distribution} angular distribution havig scale of {angular_scale}*$\pi$ for $\phi$ and {angular_scale}*$\pi/2$ for $\theta$")

    


# def basic_simulator(integrator,integator_name, experiment, engs_MeV,show_tragectories):
#     # def __one_trajectory:
        
        
        
        
        
#     if show_tragectories == True:
        
#         spec = Spectrometer_Budy(experiment.height_mm, experiment.width_mm, experiment.depth_mm,
#                                  experiment.pinhole_dia_mm, experiment.shield_mm, experiment.yoke, experiment.shield_mm)
#         fig,ax = spec._draw_spec()
    
#         # plot_field_quiver(ax,integrator)
#     experiment.solution = integrator.solution
    
#     if len(engs_MeV.shape)==1:
#         rutin = 1
#     if len(engs_MeV.shape)==2:
#         rutin = engs_MeV.shape[0]
    
#     R_vec_mm = []; v_vec_m0s = []; gamma_vec = []; label = []
    
#     for i in range(rutin):
#         # integrator.self.v0_m0s = MeV2m0s(engs_MeV[i])
#         # integrator.self.gamma = vel2gamma (integrator.self.v0_m0s)
#         if rutin == 1:
#             experiment.q_eng_MeV = engs_MeV
            
#         else:
#             experiment.q_eng_MeV = engs_MeV[i]
#         experiment._evaluate_exp_paramas()
#         params = experiment.params
        
#         integrator._update_params(params)
#         if integator_name == "Euler" or integator_name == ["Euler"]:
#             R_mm,v_m0s,gamma= integrator._euler()
#         elif integator_name == "RK2" or integator_name == ["RK2"]:
#                 R_mm,v_m0s,gamma= integrator._RK2()
#         elif integator_name == "RK4" or integator_name == ["RK4"]:
#                 R_mm,v_m0s,gamma= integrator._RK4()
#         elif integator_name == "Boris" or integator_name == ["Boris"]:
#                 R_mm,v_m0s,gamma= integrator._Boris_pusher()
                
#         R_vec_mm.append(R_mm)
#         v_vec_m0s.append(v_m0s)
#         gamma_vec.append(gamma)
        
#         if show_tragectories == True:
#             label= str(np.round(np.linalg.norm(experiment.q_eng_MeV),2))
#             pr.trajectory_plot(ax, R_mm,label)
#             print(R_mm[-1])
    
#     if show_tragectories == True:
#         # ax.legend(title='Particle total energy',title_fontsize=25,fontsize=18)
#         pos = ax.get_position()
#         ax.legend(title='Particle total energy',title_fontsize=25, loc='center left',
#           bbox_to_anchor=(1.2 , 0.5), fontsize=18)
        
#         R_vec_mm, v_vec_m0s, gamma_vec, ax
#     return R_vec_mm, v_vec_m0s, gamma_vec




    
    
    
# def real_beams(integrator,integrator_name, experiment,N_steps, total_energy_scalar_MeV,
#                number_of_experiments,energy_distribution,energy_scale, angular_distribution, angular_scale = 0.01,show_tragectories = True, show_map = True):
#     Noe = number_of_experiments
#     if energy_distribution == "Exponential":
#         total_energys = np.random.exponential(total_energy_scalar_MeV*energy_scale, Noe)
    
#     elif energy_distribution == "Gaussian":
#         total_energys = np.random.normal(loc=total_energy_scalar_MeV, scale=total_energy_scalar_MeV * energy_scale, size=Noe)
    
#     elif energy_distribution == "Uniform":
#        total_energys = np.random.uniform(low=0, high = total_energy_scalar_MeV*2, size=Noe)
      
#     elif energy_distribution == "DLA":
#          total_energys = energies_MeV = sc.stats.gamma.rvs(a=1.5, scale=total_energy_scalar_MeV, size=Noe)   
       
#     elif energy_distribution == "None":
#        total_energys = total_energy_scalar_MeV * np.ones(Noe)
    
    
#     number_of_angular_experiments = 1
#     R_engs_mat_mm =[]
#     for te in range(Noe):
#         total_energy = total_energys[te]
#         R_engs_mm = Non_collimataed_beams(integrator,integrator_name, experiment,N_steps, total_energy,number_of_angular_experiments, angular_distribution, angular_scale = angular_scale,show_tragectories = False, show_map = False)
#         R_engs_mat_mm.append(R_engs_mm)
        
#     if show_map == True:
#         fig = plt.figure()
        
#         manager = plt.get_current_fig_manager()
#         try:
#             window = manager.window
         
#             window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
#         except Exception as e:
#             print("Could not resize window:", e)
    
#         ax = fig.add_subplot(111)
#         R_engs_im_mm = np.zeros([Noe,3])
#         count = 0
#         x_exit_vec_mm = []
#         y_exit_vec_mm = []
#         for n in range(Noe):
#             R_engs_im_mm[n] = R_engs_mat_mm[n][0][-1]
#             if  abs( abs(R_engs_im_mm[n][2])-experiment.height_mm/2) < 1e-2:
#                 count +=1
#                 # x_exit_vec_mm.append(R_engs_im_mm[n,0])
#                 # y_exit_vec_mm.append(R_engs_im_mm[n,1])
#         x_exit_vec_mm = R_engs_im_mm[:,0]
#         y_exit_vec_mm = R_engs_im_mm[:,1]
#         # x_exit_vec_mm = np.array(x_exit_vec_mm
#         # y_exit_vec_mm = np.array(y_exit_vec_mm)
#         pecent = count/Noe * 100 
            
#         # if show_map==True:
            
            
#         h = plt.hist2d(x_exit_vec_mm, y_exit_vec_mm, bins=100,
#                range=[[-experiment.width_mm/2, experiment.width_mm/2], [0, experiment.depth_mm]],
#                cmap='viridis')
        
#         cbar = plt.colorbar(h[3], ax=ax)
#         cbar.set_label('counts', fontsize=25)   
#         cbar.ax.tick_params(labelsize=18)
#         plt.xlabel('Width (mm)',fontsize=25)
#         plt.ylabel('Depth (mm)',fontsize=25)
#         ax.tick_params(axis='both', labelsize=18)
#         ax.axis('equal')
#         ax.set_title(f"{pecent} of the particle passed the ditector with {angular_distribution} angular distribution havig scale of {angular_scale}*$\pi$ for $\phi$ and {angular_scale}*$\pi/2$ for $\theta$")
    
#     R_engs_return_mm = []
#     count = 0
#     for n in range(Noe):
#         R_engs_return_mm.append(R_engs_mat_mm[n][0][-1])
#         if  abs( abs(R_engs_return_mm[n][2])-experiment.height_mm/2) < 1e-2:
#             count +=1
#             # x_exit_vec_mm.append(R_engs_im_mm[n,0])
#             # y_exit_vec_mm.append(R_engs_im_mm[n,1])
#     pecent = count/Noe * 100 
#     return np.array(R_engs_return_mm), pecent