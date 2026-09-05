#from IPython import get_ipython
#import matplotlib.pyplot as plt
#plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');


import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

from spectrometer.geometry import geometry
import spectrometer.integrators as integ
import spectrometer.plot_result as pr
import spectrometer.fields as fields




if __name__ == "__main__":
    me_kg = 9.109*1e-31
    B_T = np.array([0.5,0.0,0.0])
    E_V0m = np.array([0,0,0])
    R0_mm  = np.array([0,-12.5,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    steps = 150
    yoke = 1
    shield_mm = 12.5
    ig,ax  = geometry(shield_mm = shield_mm, yoke=yoke)
    pinhole_dia_mm = 3
    fringe = 1
    Z_mm = integ.analitic_sol_vel2dist (e_eng_MeV[1], me_kg, e_C, height_mm, B_T[0])
    print(Z_mm)
    
    R_vec_mm,v_vec_m0s,gamma_vec = integ.euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    pr.trajectory_plot(ax, R_vec_mm)
    print(R_vec_mm[-1])
    
    R_vec_mm,v_vec_m0s,gamma_vec = integ.RK2 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    pr.trajectory_plot(ax, R_vec_mm)
    print(R_vec_mm[-1])
    
    R_vec_mm,v_vec_m0s,gamma_vec = integ.RK4 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    pr.trajectory_plot(ax, R_vec_mm)
    print(R_vec_mm[-1])
    
    R_vec_mm,v_vec_m0s,gamma_vec = integ.Boris_pusher (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    pr.trajectory_plot(ax, R_vec_mm)
    print(R_vec_mm[-1])
    
    
    mag = fields.Magenetic_field_Analitic(B_T[0], width_mm, depth_mm, k=0, fringe=fringe, sharp_edge=1, yoke=yoke)
    
    fields.plot_magnetic_quiver(ax, mag, width_mm, depth_mm, height_mm, yoke,fringe, shield_mm)
    