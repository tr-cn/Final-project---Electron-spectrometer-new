import numpy as np
me_kg = 9.109*1e-31
B_T = np.array([0,0.5,0])
E_V0m = 0
R0 = np.array([0,0,0])
e_C = 1.602*1e-19
e_eng_MeV = 10
height_mm =26; width_mm = 12.5; depth_mm   = 50.8


def analitic_sol_vel2dist (q_eng_MeV, m_kg, q_C, height_mm, B_T):
    # An explanation of how radius and the velocities are calculated is given
    # in the documation
    h_mm = height_mm/2
    q_eng_J = q_eng_MeV*1e6 * 1.602*1e-19
    c_m0s = 299792458
    v_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
    gamma = 1 / np.sqrt ( 1 - (v_m0s/c_m0s)**2 )
    R_m = gamma * m_kg*v_m0s / (abs(q_C) * B_T)
    R_mm = R_m *1e3
    
    Z_mm = np.sqrt ( 2*R_mm*h_mm - h_mm**2 )
    return Z_mm

#z = analitic_sol_vel2dist (e_eng_MeV, me_kg, e_C, height_mm , B_T)


def vel2gamma (v_m0s):
    c_m0s = 299792458
    v0_mag_m0s = np.linalg.norm(v_m0s)
    return 1 / np.sqrt ( 1 - (v0_mag_m0s/c_m0s)**2 )

def euler (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0):
    h_mm = height_mm/2; d_mm = depth_mm; w_mm = width_mm/2;
    
    q_eng_J = q_eng_MeV*1e6 * 1.602*1e-19
    c_m0s = 299792458
    v0_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
    gamma = vel2gamma(v0_m0s);
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_mm = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_mm.append(R0);
    
    CFL = 0.5; dx_m = 1e-6; dt_s = dx_m*CFL; # not realy neaded in euler
    while R_vec_mm[-1][2] < h_mm and R_vec_mm[-1][1] < d_mm and abs(R_vec_mm[-1][0]) < abs (w_mm):
        
        v_m0s = np.array(v_vec_m0s[-1]); 
        dv_m0s = E_V0m*dt_s + q_C*np.cross(v_m0s,B_T)
        dr_m = dv_m0s*dt_s;
        
        v_new_m0s = v_m0s + dv_m0s;   gamma = vel2gamma(v_m0s);
        v_fixed_m0s = v_new_m0s*gamma
        
        v_vec_m0s.append(v_fixed_m0s)
        R_vec_mm.append(R_vec_mm[-1] + dr_m*1e-3)
        gamma_vec.append(gamma)
        
        
    return R_vec_mm,v_vec_m0s,gamma_vec
        
        
        


R_vec_mm,v_vec_m0s,gamma_vec = euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0)

    
    
    
    

