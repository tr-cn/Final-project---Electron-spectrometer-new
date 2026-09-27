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
    
    
    
    
    # e_engs_y_MeV = [i for i in np.linspace(0.5, 15,4)]
    # e_engs_y_MeV = 0.6

    # e_engs_MeV = np.array([ [0,i,0] for i in np.linspace(0.6, 15,4)])
    # e_engs_MeV = np.array([0,0.3,0])
    
    # integrator_name = "Boris"#"Euler"#"RK4"#"RK2"
    
    # basic_simulator(integrator,integrator_name, experiment, e_engs_MeV,show_tragectories = True)
    
    
    integrator_name = "Boris" #"Boris"#"Euler"#"RK4"#"RK2"
    e_engs_MeV = np.array([ [0,i,0] for i in np.linspace(0.6, 15,4)])
    N_steps_vec = [5e3,1e4,5e4]
    # basic_simulator(integrator,integrator_name, experiment, e_engs_MeV,show_tragectories = True)
    integrators_name = ["Boris","RK2","RK4","Euler"]
    energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
    
    
    
    
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
    