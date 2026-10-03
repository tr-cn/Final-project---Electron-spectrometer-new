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
    R0_mm  = np.array([0,-15,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = -5e7
    yoke = 0
    shield_mm = 5
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.05
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =    "Analitic field" #  "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**8 + 1;
    Ny_p = 2**8 + 1;
    

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
    
    # %%    #_________________Integration environment intitiation__________________
    if True:#False:
        params["solution"] = "Analitic field"
    
        integrator = Integrators(params)
    
    if False:#True:
        params["solution"] = "Numeric Field"
        integrator = Integrators(params)
        # integrator.psi._show_potential()
        # integrator.mag._show_field()
        # integrator.phi._show_potential()
        # integrator.Ele._show_field()




# %% # _______________________TNSA - experiment_________________________
    if True:
        integrator_name = ["RK4_Herm"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 1
        params["Ex0_Vm"] = 5e6
        params["R0_mm"] = np.array([0,-5,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        
        
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential"#"DLA"
        angular_scale = 0.005
        energy_scale = 0.1
        e_engs_MeV_scalars = [3,0.1]
        mq_vec_kg = [9.109*1e-31,1.67*1e-27]
        q_vec_C = [-1.602*1e-19,-1.602*1e-19]
 
        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalars_MeV = e_engs_MeV_scalars, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
        


N_steps_vec = None#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential"#"DLA"
        angular_scale = 0.005
        energy_scale = 0.1
        e_engs_MeV_scalars = [2,2]
        mq_vec_kg  =[9.109*1e-31, 5*9.1*1e-31]# [9.1e-30]#[9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,0.5*-1.602*1e-19]#[-1.602*1e-18]#[-1.602*1e-19,-1.602*1e-18]
 
        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalars_MeV = e_engs_MeV_scalars, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)



    if True:
        params = copy.deepcopy(params_beckup)
        params = experiment.params 
        integrator_name = ["RK4_Herm"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 1e7
        params["R0_mm"] = np.array([0,0,0])
        params["N_steps"] = 5e2
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        
        
        N_steps_vec = None#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential"#"DLA"
        angular_scale = 0.005
        energy_scale = 0.1
        e_engs_MeV_scalars = [2,2]
        mq_vec_kg  =[9.109*1e-31, 5*9.1*1e-31]# [9.1e-30]#[9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,0.5*-1.602*1e-19]#[-1.602*1e-18]#[-1.602*1e-19,-1.602*1e-18]
 
        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalars_MeV = e_engs_MeV_scalars, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)