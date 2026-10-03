#from IPython import get_ipython
#import matplotlib.pyplot as plt
#plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');


# %%
import matplotlib.pyplot as plt

import numpy as np
from types import SimpleNamespace

plt.close('all')

from spectrometer.geometry import Spectrometer_Budy
from spectrometer.Experiment import Experiment
from spectrometer.integrators import Integrators
# import spectrometer.integrators as integ
import spectrometer.plot_result as pr
import spectrometer.analitic_fields as analitic_fields
from spectrometer.simulators import *
import copy




if __name__ == "__main__":
   
    
    me_kg =9.109*1e-31 #9.109*1e-31
    R0_mm  = np.array([0,-10,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = 1e8

    yoke = 0
    shield_mm = 5
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.05
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =    "Numeric Field"  #  "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**4 + 1;
    Ny_p = 2**4 + 1;
    
    grid_interpulation = "Linear" # "Spline"#
    # spec = Spectrometer_Budy(shield_mm = shield_mm, yoke=yoke)
    # fig,ax = spec._draw_spec()
    
    experiment = Experiment(
                            height_mm = height_mm, width_mm = width_mm, depth_mm = depth_mm, shield_mm = shield_mm, pinhole_dia_mm = pinhole_dia_mm, yoke = yoke, # Spectrometere and simulation border
                            R0_mm= R0_mm, q_eng_MeV = e_eng_MeV, m_kg = me_kg, q_C = e_C, # Particle parameters
                            Bx0_T = Bx0_T, Ex0_Vm = Ex0_Vm, fringe = fringe, sharp_edge = sharp_edge, # Fildes parameters
                            CFL = CFL, N_steps = N_steps, solution = solution, # Simulation resolutions
                            Nx_p = Nx_p, Ny_p = Ny_p, grid_interpulation = grid_interpulation, # Potential solving
                            )
       
    experiment._evaluate_exp_paramas()
    params = experiment.params
   
    # %%    #_________________Integration environment intitiation__________________
    if False:#False:
        params["solution"] = "Analitic field"
    
        integrator = Integrators(params)
    
    if True:#True:
        params["solution"] = "Numeric Field"
        integrator = Integrators(params)
        
        # integrator.psi._show_potential()
        # integrator.mag._show_field()
        # integrator.phi._show_potential()
        # integrator.Ele._show_field()
        
    params_beckup = copy.deepcopy(params)
    #%%    Fully analitic assuming no velocity in Z direction, no fringe fields
    if False:
        params = copy.deepcopy(params_beckup)
        params["Ex0_Vm"] = 5e7
        params["R0_mm"] = np.array([0,-4,0])
        
        params = experiment._update_params(params)
        # integrator.mag.interpolation
        params = integrator._update_params(params)
        e_engs_vec_MeV = np.array([ [0,i,0] for i in np.linspace(1,15,1001)])
        m_vec_kg  =[9.109*1e-31, 5*9.1*1e-31]# [9.1e-30]#[9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,0.5*-1.602*1e-19]#[-1.602*1e-18]#[-1.602*1e-19,-1.602*1e-18]
        # params["q_eng_MeV"] = np.array([10,10,0])
        Fully_analitic(integrator,experiment,e_engs_vec_MeV,m_vec_kg,q_vec_C)
        
        
        
        
        
    # %%    #___________________________Basic simuations_______________________
    
    
    
    if True:
        
        e_engs_MeV = np.array([ [0,i,0] for i in np.linspace(1,15,15)])
        # e_engs_MeV = np.array([-0.005,2,0.005])# When using this, it is better to work with velocity in the y driection so the energy will not exceed the speed of light
        
        integrator_name = "Boris"#"Euler"#"RK4"#"RK2"
        basic_simulator(integrator,integrator_name, experiment, e_engs_MeV,show_tragectories = True)
    

    # %%    #___________________________Converssion Tests______________________
    if True:
        params = copy.deepcopy(params_beckup)
        integrators_name =["Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]#["RK4_Herm"]#
        params["fringe"] = 1
        params["Ex0_Vm"] = Ex0_Vm
        params["R0_mm"] = np.array([0,-4,0])
        params["grid_interpulation"] = "Spline"#"Spline"# "Linear"
        params = experiment._update_params(params)
        # integrator.mag.interpolation
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        e_engs_MeV = np.array([0, 3, 0])
        # N_steps_vec = [1e2,5e2,1e3,1e4]
# 
        N_steps_vec = [25,50,100,200,400,800,1600,3200,6400]

        convergence_test (integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
 # %%    # _____________Energy Vs N_steps Vs Integrator_____________________
    if False:
        params = copy.deepcopy(params_beckup)
        integrators_name =["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 0
        params["R0_mm"] = np.array([0,0,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        e_engs_MeV = np.array([0,5,0])
        N_steps_vec = [1e2]#,4e3,5e3,1e4,3e4,5e4]
         
        energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
        
    # %%# _____________________Non_collimataed_beams___________________________
    
    if False:
        params = copy.deepcopy(params_beckup)
        integrators_name = ["RK4_Herm"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 1
        params["Ex0_Vm"] = 0
        params["R0_mm"] = np.array([0,-3,0])
        params["N_steps"] = 5e2
        params = experiment._update_params(params)
        integrator._update_params(params)
        # integrator._update_solution(integrator.solution)
        
        e_engs_MeV_scalar = 4
        N_steps_vec = None#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Exponential"#"Uniform" #"Gaussian" #"Exponential"
        Non_collimataed_beams(integrator,integrators_name, experiment,N_steps_vec, e_engs_MeV_scalar, number_of_experiments, angular_distribution,angular_scale=0.03, show_tragectories = False, show_map = True)
        
    # %% # ______________Real_beams - spectral distribution__________________ 
    if False:
        
        integrator_name = ["RK4_Lin"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 1e6
        params["R0_mm"] = np.array([0,-5,0])
        params["N_steps"] = 5e2
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        e_engs_MeV_scalar = 4
        N_steps_vec = None#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential" 
        angular_scale = 0.005
        energy_scale = 0.1
        real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
                   N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
    # %% # _______________________TNSA - experiment_________________________
    if True:
        params = copy.deepcopy(params_beckup)
        params = experiment.params 
        integrator_name = ["RK4_Herm"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 5e7#1e7
        params["R0_mm"] = np.array([0,0,0])
        params["N_steps"] = 5e2
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        
        
        N_steps_vec = None#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 100
        angular_distribution = "None"# "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "Uniform"#"DLA" #"None"#"Uniform" #"Gaussian" #"Exponential"#"DLA"
        angular_scale = 0.005
        energy_scale = 0.1
        e_engs_MeV_scalars = [15,7]#[2,2]
        mq_vec_kg  =[9.109*1e-31, 5*9.1*1e-31]# [9.1e-30]#[9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,0.5*-1.602*1e-19]#[-1.602*1e-18]#[-1.602*1e-19,-1.602*1e-18]
 
        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalars_MeV = e_engs_MeV_scalars, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
        
        
            
#%%        
    # if True:
    #     integrator_name = ["Boris"]#["Boris","RK2","RK4","Euler"]
    #     e_engs_MeV_scalar = 4
    #     N_steps_vec = 1e3#,4e3,5e3,1e4,3e4,5e4]
    #     number_of_experiments = 60
    #     angular_distribution = "Uniform" #"Gaussian" #"Exponential"
    #     energy_distribution =  "Uniform" #"Gaussian" #"Exponential"
    #     angular_scale = 0.03
    #     energy_scale = 0.1
    #     real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
    #                N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
    #                number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
    #                angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
        
        
        
# def dNdE_shape(E, Te):
#     return np.sqrt(E/np.pi) * Te**(-1.5) * np.exp(-E/Te)

# Te = 3.0   # טמפרטורת-האלקטרונים, MeV
# E_peak = Te / 2                              # שיא הפונקציה (נגזרת אנליטית)
# f_max = dNdE_shape(E_peak, Te)
# E_max_range = 15 * Te                        # תחום גזירה, כולל ה"זנב" האקספוננציאלי

# def rejection_sample_energy(N, Te, E_max_range, f_max):
#     samples = []
#     while len(samples) < N:
#         E_candidates = np.random.uniform(0, E_max_range, N)
#         y_candidates = np.random.uniform(0, f_max, N)
#         accepted = E_candidates[y_candidates <= dNdE_shape(E_candidates, Te)]
#         samples.extend(accepted.tolist())
#     return np.array(samples[:N])

# energies_MeV = rejection_sample_energy(N=1000, Te=Te, E_max_range=E_max_range, f_max=f_max)
    
    
    
    # simulation_engs_collimated(integrator_anal, experiment, e_engs_MeV,ax)
    
    
    
        
    
    
    # integrator.solution = "Numeric Field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._RK4 ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # integrator.solution = "Analitic field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._RK4 ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    
    
    # integrator.solution = "Numeric Field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._Boris_pusher ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # integrator.solution = "Analitic field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._Boris_pusher ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    
    
    
    
    
    
    
    # Z_mm = integ.analitic_sol_vel2dist (e_eng_MeV[1], me_kg, e_C, height_mm, B_T[0])
    # print(Z_mm)
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.RK2 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.RK4 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.Boris_pusher (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    # mag = fields.Magenetic_field_Analitic(B_T[0], width_mm, depth_mm, k=0, fringe=fringe, sharp_edge=1, yoke=yoke)
    
    # fields.plot_magnetic_quiver(ax, mag, width_mm, depth_mm, height_mm, yoke,fringe, shield_mm)
    