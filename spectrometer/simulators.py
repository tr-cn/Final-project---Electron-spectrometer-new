import numpy as np
import matplotlib.pyplot as plt
from spectrometer.Experiment import MeV2m0s, vel2gamma
from spectrometer.Experiment import Experiment
import spectrometer.plot_result as pr
from spectrometer.geometry import *


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
            eps = 10**abs(round(np.log10(integrator.Bx0_T),0))
        if r == 1:
            eps = 10**abs(round(np.log10(integrator.Ex0_Vm),0))
        Fx, Fy, Fz = _log_scale_vectors(Fx, Fy, Fz, eps=eps)
        ax.quiver(X_mm, Y_mm, Z_mm, 
                  Fx, Fy, Fz, 
                  length=quive_len,       
                  normalize = False,
                  color=color, 
                  alpha=0.6, 
                  arrow_length_ratio = 1)
    
    
    
def basic_simulator(integrator,integator_name, experiment, engs_MeV,show_tragectories):

    if show_tragectories == True:
        spec = Spectrometer_Budy(experiment.height_mm, experiment.width_mm, experiment.depth_mm,
                                 experiment.pinhole_dia_mm, experiment.shield_mm, experiment.yoke, experiment.shield_mm)
        fig,ax = spec._draw_spec()
    
        # plot_field_quiver(ax,integrator)
    experiment.solution = integrator.solution
    
    if len(engs_MeV.shape)==1:
        rutin = 1
    if len(engs_MeV.shape)==2:
        rutin = engs_MeV.shape[0]
    
    R_vec_mm = []; v_vec_m0s = []; gamma_vec = [];
    
    for i in range(rutin):
        # integrator.self.v0_m0s = MeV2m0s(engs_MeV[i])
        # integrator.self.gamma = vel2gamma (integrator.self.v0_m0s)
        if rutin == 1:
            experiment.q_eng_MeV = engs_MeV
        else:
            experiment.q_eng_MeV = engs_MeV[i]
        experiment._evaluate_exp_paramas()
        params = experiment.params
        
        integrator._update_params(params)
        if integator_name == "Euler":
            R_mm,v_m0s,gamma= integrator._euler()
        elif integator_name == "RK2":
                R_mm,v_m0s,gamma= integrator._RK2()
        elif integator_name == "RK4":
                R_mm,v_m0s,gamma= integrator._RK4()
        elif integator_name == "Boris":
                R_mm,v_m0s,gamma= integrator._Boris_pusher()
                
        R_vec_mm.append(R_mm)
        v_vec_m0s.append(v_m0s)
        gamma_vec.append(gamma)
        
        if show_tragectories == True:
            pr.trajectory_plot(ax, R_mm)
            print(R_mm[-1])
    if show_tragectories == True:
        R_vec_mm, v_vec_m0s, gamma_vec, ax
    return R_vec_mm, v_vec_m0s, gamma_vec





    


def energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, engs_MeV,N_steps_vec):
    
    params = experiment.params
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    colors = ["magenta", "red", "blue", "green"]
    for i in range (len(integrators_name)):
        integrator_name = integrators_name[i]
        R_vec_mm = [];  v_vec_m0s = []; gamma_vec = []; gamma_retio_vec = [] 
        
        for n in range(len(N_steps_vec)):
            # print (N_steps_vec[n])
            params["N_steps"] = N_steps_vec[n]
            params = experiment._update_params(params)
            integrator._update_params(params)
            R_engs_mm, v_engs_m0s, gamma_engs_vec = basic_simulator (integrator,integrator_name, experiment, engs_MeV,show_tragectories=False)
            gamma_retio = []
            norm_engs_MeV = []
            for g in range(len(gamma_engs_vec)):
                # print(gamma_engs_vec[n])
                norm_engs_MeV.append(np.linalg.norm(engs_MeV[g]))
                gamma_retio.append(gamma_engs_vec[g][-1]/gamma_engs_vec[g][0]*norm_engs_MeV[-1])
            
            
            R_vec_mm.append(R_engs_mm);
            v_vec_m0s.append(v_engs_m0s);
            gamma_vec.append(gamma_engs_vec)
            gamma_retio_vec.append(gamma_retio)
            norm_engs_MeV = np.array(norm_engs_MeV)
            
        gamma_MeV = np.array(gamma_retio_vec)
        
         
        # color = colors[i]
        # ax.plot(N_steps_vec, gamma_MeV,color=color)    
        color = colors[i]
        for g in range(gamma_MeV.shape[1]):
            
            label = integrator_name if g == 0 else None   # legend entry אחד בלבד לכל אינטגרטור
            ax.scatter(N_steps_vec, gamma_MeV[:, g], color=color, alpha=0.6, s=50,
                       edgecolor='black', linewidth=0.5, label=label)

    ax.set_xlabel('N_steps')
    ax.set_ylabel('gamma ratio')
    ax.set_xscale('log')
    ax.legend(title='Integrator')
    
        
        
        
        
        
        
        