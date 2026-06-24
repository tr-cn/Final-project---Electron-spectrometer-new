#from IPython import get_ipython
#import matplotlib.pyplot as plt
#plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');


import matplotlib.pyplot as plt
import numpy as np


from spectrometer.geometry import geometry
import spectrometer.physics as phys
import spectrometer.plot_result as pr





if __name__ == "__main__":
    me_kg = 9.109*1e-31
    B_T = np.array([0.5,0,0])
    E_V0m = np.array([0,0,0])
    R0 = np.array([0,0,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8


    ig,ax  = geometry()
    R_vec_mm,v_vec_m0s,gamma_vec = phys.euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0)
    
    pr.trajectory_plot(ax, R_vec_mm)