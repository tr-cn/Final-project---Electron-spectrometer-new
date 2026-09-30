#from IPython import get_ipython
#import matplotlib.pyplot as plt
#plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');


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




if __name__ == "__main__":
   
    
    me_kg = 9.109*1e-31
    R0_mm  = np.array([0,0,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = 0*1e8
    yoke = 1
    shield_mm = 0
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.1
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =   "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**5 + 1;
    Ny_p = 2**5 + 1;
    

    # spec = Spectrometer_Budy(shield_mm = shield_mm, yoke=yoke)
    # fig,ax = spec._draw_spec()
    
    experiment = Experiment(
                            height_mm = height_mm, width_mm = width_mm, depth_mm = depth_mm, shield_mm = shield_mm, pinhole_dia_mm = pinhole_dia_mm, yoke = yoke, # Spectrometere and simulation border
                            R0_mm= R0_mm, q_eng_MeV = e_eng_MeV, m_kg = me_kg, q_C = e_C, # Particle parameters
                            Bx0_T = Bx0_T, Ex0_Vm = Ex0_Vm, fringe = fringe, sharp_edge = sharp_edge, # Fildes parameters
                            CFL = CFL, N_steps = N_steps, solution = solution, # Simulation resolutions
                            Nx_p = Nx_p, Ny_p = Ny_p # Potential solving
                            )
       
    experiment._evaluate_exp_paramas()
    params = experiment.params
    

    if True:#False:
        params["solution"] = "Analitic field"
    
        integrator = Integrators(params)
    
    if False:#True:
        params["solution"] = "Numeric Field"
        integrator = Integrators(params)
    
    
    
    
    if False:
        e_engs_MeV = np.array([ [0,i,0] for i in np.linspace(10,20,4)])
        # e_engs_MeV = np.array([-0.005,2,0.005])# When using this, it is better to work with velocity in the y driection so the energy will not exceed the speed of light
        
        integrator_name = "Boris"#"Euler"#"RK4"#"RK2"
        basic_simulator(integrator,integrator_name, experiment, e_engs_MeV,show_tragectories = True)
    

   ###########################################################################
    # __________Energy Vs N_steps Vs Integrator__________________
    if False:
        integrators_name = ["RK4","Boris","RK2"]#["Boris","RK2","RK4","Euler"]
        e_engs_MeV = np.array([[0,6,0],[0,10,0],[0,16,0]])
        N_steps_vec = [1e3,2e3,3e3]#,4e3,5e3,1e4,3e4,5e4]
        
        energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
        
###############################################################################
    
    if False:
        integrators_name = ["Boris"]#["Boris","RK2","RK4","Euler"]
        e_engs_MeV_scalar = 10
        N_steps_vec = 1e3#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1
        angular_distribution = "Exponential"#"Uniform" #"Gaussian" #"Exponential"
        Non_collimataed_beams(integrator,integrators_name, experiment,N_steps_vec, e_engs_MeV_scalar, number_of_experiments, angular_distribution,angular_scale=0.03, show_tragectories = False, show_map = True)
        
        
    if False:
        integrator_name = ["Boris"]#["Boris","RK2","RK4","Euler"]
        e_engs_MeV_scalar = 4
        N_steps_vec = 1e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 60
        angular_distribution = "None"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential" 
        angular_scale = 0.03
        energy_scale = 0.1
        real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
                   N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
    if True:
        Ex0_Vm = 1e7
        params["Ex0_Vm"] = Ex0_Vm
        experiment._update_params(params)
        integrator._update_params(params)
        
        
        
        integrator_name = ["Boris"]#["Boris","RK2","RK4","Euler"]
        e_engs_MeV_scalar = 4
        N_steps_vec = 1e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 3
        angular_distribution = "None"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential" 
        angular_scale = 0.03
        energy_scale = 0.1
        
        mq_vec_kg = [9.109*1e-31,9.109*1e-30]
        q_vec_C = [-1.602*1e-19,-1.602*1e-16]


        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
        
        
            
        
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
    