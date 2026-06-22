import numpy as np
me_kg = 9.109*10**-31
B_T = 0.5
e_C = 1.602*10**-19
e_eng_MeV = 10
height_mm = 26

def analitic_sol_vel2dist (q_eng_MeV, m_kg, q_C, height_mm, B_T):
    # An explanation of how radius and the velocities are calculated is given
    # in the documation
    h_mm = height_mm/2
    c_m0s = 299792458
    v_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_MeV + m_kg*c_m0s**2) )**2 )# [m/s]
    R_m = m_kg*v_m0s / (abs(q_C) * B_T)
    R_mm = R_m *1**-3
    
    Z_mm = np.sqrt (2*R_mm*h_mm - h_mm**2 )
    
    return Z_mm


analitic_sol_vel2dist (e_eng_MeV, me_kg, e_C, height_mm , B_T)