import numpy as np
from spectrometer.Experiment import MeV2m0s, vel2gamma
from spectrometer.Experiment import Experiment
import spectrometer.plot_result as pr

def simulation_engs_collimated(integrator,experiment, engs_MeV,ax):
    # spec = Spectrometer_Budy(shield_mm = shield_mm, yoke=yoke)
    # fig,ax = spec._draw_spec()
    
    for i in range(len(engs_MeV)):
        # integrator.self.v0_m0s = MeV2m0s(engs_MeV[i])
        # integrator.self.gamma = vel2gamma (integrator.self.v0_m0s)
        experiment.q_eng_MeV = engs_MeV[i]
        experiment._evaluate_exp_paramas()
        params = experiment.params
        
        integrator._update_params(params)
        R_vec_mm,v_vec_m0s,gamma_vec = integrator._RK2()
        pr.trajectory_plot(ax, R_vec_mm)
        print(R_vec_mm[-1])
        # print(integrator_num.q_eng_MeV)
