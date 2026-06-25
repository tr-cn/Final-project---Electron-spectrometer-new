import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#plt.rcParams.update({
#   'font.family': 'serif',
#  'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
# 'mathtext.fontset': 'stix'  # Matches Times New Roman math styles seamlessly
#})

from matplotlib.ticker import MaxNLocator, FormatStrFormatter






def analitic_sol_vel2dist (q_eng_MeV, m_kg, q_C, height_mm, B_T):
    # An explanation of how radius and the velocities are calculated is given
    # in the documation
    q_C = abs(q_C)
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


def MeV2m0s(q_eng_MeV, m_kg):
    q_eng_J = abs(q_eng_MeV)*1e6 * 1.602*1e-19
    c_m0s = 299792458
    v0_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
    return v0_m0s

def vel2gamma (v_m0s):
    c_m0s = 299792458
    v0_mag_m0s = np.linalg.norm(v_m0s)
    return 1 / np.sqrt ( 1 - (v0_mag_m0s/c_m0s)**2 )



def get_magnetic_field (B_type, dx_mm, dy_mm, dz_mm): # Need to think how to do it correctly
   
    return


def get_electric_field (E_type, dx_mm, dy_mm, dz_mm): # Need to think how to do it correctly
   
    return


def euler (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0):
    c_m0s = 299792458
    h_mm = height_mm/2; d_mm = depth_mm; w_mm = width_mm/2;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg)
    
    gamma = vel2gamma(v0_m0s)
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_mm = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_mm.append(R0);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/10000
    while R_vec_mm[-1][2] < h_mm and R_vec_mm[-1][1] < d_mm and R_vec_mm[-1][1]>=0 and abs(R_vec_mm[-1][0]) < abs (w_mm):
        gamma = vel2gamma(v_vec_m0s[-1])
        v_m0s = np.array(v_vec_m0s[-1]); 
        dv_dt_m0s2= (q_C/(gamma*m_kg)) * (E_V0m + np.cross(v_m0s,B_T) -  v_m0s*np.dot(E_V0m, v_m0s) /c_m0s**2 )
        
        v_new_m0s = v_m0s + dv_dt_m0s2*dt_s
        
        dr_m = v_new_m0s*dt_s;
        
        v_vec_m0s.append(v_new_m0s)
        R_vec_mm.append(R_vec_mm[-1] + dr_m*1e3)
        gamma_vec.append(vel2gamma(v_new_m0s))
        
    return R_vec_mm,v_vec_m0s,gamma_vec
        


#R_vec_mm,v_vec_m0s,gamma_vec = euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0)


    

def RK2 (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0):
    c_m0s = 299792458
    h_mm = height_mm/2; d_mm = depth_mm; w_mm = width_mm/2;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg)
    
    gamma = vel2gamma(v0_m0s)
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_mm = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_mm.append(R0);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/10000
    while R_vec_mm[-1][2] < h_mm and R_vec_mm[-1][1] < d_mm and R_vec_mm[-1][1]>=0 and abs(R_vec_mm[-1][0]) < abs (w_mm):
        v_m0s_i = np.array(v_vec_m0s[-1]); 
        gamma_i = vel2gamma(v_vec_m0s[-1])
        R_m_i = np.array(R_vec_mm[-1])/1e3; 
       
        dv_dt_m0s2_i = (q_C/(gamma_i*m_kg)) * (E_V0m + np.cross(v_m0s_i,B_T) -  v_m0s_i*np.dot(E_V0m, v_m0s_i) /c_m0s**2 )
        dr_dt_m0s_i = v_m0s_i;
        
        k1_v = 1/2*dv_dt_m0s2_i*dt_s
        k1_r = 1/2*dr_dt_m0s_i*dt_s
        
        
        v_m0s_m = v_m0s_i + k1_v; 
        gamma_m = vel2gamma(v_m0s_m)
        R_m_m = R_m_i + k1_r; # not realy neaded
        
        dv_dt_m0s2_m = (q_C/(gamma_m*m_kg)) * (E_V0m + np.cross(v_m0s_m,B_T) -  v_m0s_m*np.dot(E_V0m, v_m0s_m) /c_m0s**2 )
        dr_dt_m0s_m = v_m0s_m;
        
        
        k2_v = dt_s * dv_dt_m0s2_m
        k2_r = dt_s * dr_dt_m0s_m
        
        
        v_new_m0s = v_m0s_i + k2_v
        r_new_mm = (R_m_i + k2_r)*1e3
        gamma_new = vel2gamma(v_new_m0s)
        
        
        v_vec_m0s.append(v_new_m0s)
        R_vec_mm.append(r_new_mm)
        gamma_vec.append(gamma_new)
        
    return R_vec_mm,v_vec_m0s,gamma_vec
        




def RK4 (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0):
    c_m0s = 299792458
    h_mm = height_mm/2; d_mm = depth_mm; w_mm = width_mm/2;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg)
    
    gamma = vel2gamma(v0_m0s)
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_mm = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_mm.append(R0);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/10000
    while R_vec_mm[-1][2] < h_mm and R_vec_mm[-1][1] < d_mm and R_vec_mm[-1][1]>=0 and abs(R_vec_mm[-1][0]) < abs (w_mm):
        v_m0s_i = np.array(v_vec_m0s[-1]); 
        gamma_i = vel2gamma(v_vec_m0s[-1])
        R_m_i = np.array(R_vec_mm[-1])/1e3; 
       
        dv_dt_m0s2_i = (q_C/(gamma_i*m_kg)) * (E_V0m + np.cross(v_m0s_i,B_T) -  v_m0s_i*np.dot(E_V0m, v_m0s_i) /c_m0s**2 )
        dr_dt_m0s_i = v_m0s_i;
        
        k1_v = dv_dt_m0s2_i*dt_s
        k1_r = dr_dt_m0s_i*dt_s
        
                
        v_m0s_m1 = v_m0s_i + k1_v/2; 
        gamma_m1= vel2gamma(v_m0s_m1)
        R_m_m1 = R_m_i + k1_r/2; # not realy neaded
        
        
        dv_dt_m0s2_m1 = (q_C/(gamma_m1*m_kg)) * (E_V0m + np.cross(v_m0s_m1,B_T) -  v_m0s_m1*np.dot(E_V0m, v_m0s_m1) /c_m0s**2 )
        dr_dt_m0s_m1 = v_m0s_m1;        
             
        
        k2_v = dv_dt_m0s2_m1*dt_s
        k2_r = dr_dt_m0s_m1*dt_s
        
        v_m0s_m2 = v_m0s_i + k2_v/2; 
        gamma_m2= vel2gamma(v_m0s_m2)
        R_m_m2 = R_m_i + k2_r/2; # not realy neaded
        
        dv_dt_m0s2_m2 = (q_C/(gamma_m2*m_kg)) * (E_V0m + np.cross(v_m0s_m2,B_T) -  v_m0s_m2*np.dot(E_V0m, v_m0s_m2) /c_m0s**2 )
        dr_dt_m0s_m2 = v_m0s_m2;   
        
        
        k3_v = dv_dt_m0s2_m2*dt_s
        k3_r = dr_dt_m0s_m2*dt_s
        
        
        v_m0s_m3 = v_m0s_i + k3_v; 
        gamma_m3= vel2gamma(v_m0s_m3)
        R_m_m3 = R_m_i + k3_r; # not realy neaded
        
        dv_dt_m0s2_m3 = (q_C/(gamma_m3*m_kg)) * (E_V0m + np.cross(v_m0s_m3,B_T) -  v_m0s_m3*np.dot(E_V0m, v_m0s_m3) /c_m0s**2 )
        dr_dt_m0s_m3 = v_m0s_m3;  
        
        k4_v = dv_dt_m0s2_m3*dt_s
        k4_r = dr_dt_m0s_m3*dt_s 
        
        
        v_new_m0s = v_m0s_i + 1/6 * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        r_new_mm = (R_m_i + 1/6 * (k1_r + 2*k2_r + 2*k3_r + k4_r))*1e3
        gamma = vel2gamma(v_new_m0s)
        
        
        v_vec_m0s.append(v_new_m0s)
        R_vec_mm.append(r_new_mm)
        gamma_vec.append(gamma)
        
    return R_vec_mm,v_vec_m0s,gamma_vec
        




def Boris_pusher (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0):
    c_m0s = 299792458
    h_mm = height_mm/2; d_mm = depth_mm; w_mm = width_mm/2;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg)
    
    gamma = vel2gamma(v0_m0s)
    
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/10000
    
    v_minus_half_m0s = v0_m0s - 1/2 *( (q_C/(gamma*m_kg)) * (E_V0m + np.cross(v0_m0s,B_T) - v0_m0s*np.dot(E_V0m, v0_m0s) /c_m0s**2 ) )* dt_s
    
    v_minus_half_vec_m0s = list([]); gamma_vec = list([]); R_vec_mm = list([]); 
    v_minus_half_vec_m0s.append(v_minus_half_m0s); gamma_vec.append(gamma); R_vec_mm.append(R0);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    
    
       
    
    while abs(R_vec_mm[-1][2])  < abs(h_mm) and R_vec_mm[-1][1] < d_mm and R_vec_mm[-1][1]>=0 and abs(R_vec_mm[-1][0]) < abs (w_mm):
        
        v_minus_half_m0s = np.array(v_minus_half_vec_m0s[-1]); 
        gamma_i = vel2gamma(v_minus_half_m0s)
        R_m_i = np.array(R_vec_mm[-1])/1e3; 
        
        
        v1 = v_minus_half_m0s + (q_C/(gamma_i*m_kg))*  (E_V0m  - v_minus_half_m0s*np.dot(E_V0m, v_minus_half_m0s) /c_m0s**2) * dt_s/2
        
        gamma_m1 =  vel2gamma(v1)
        v2 = np.zeros_like(v1);
        v2[0] = v1[0];
        wc = abs(q_C*np.linalg.norm(B_T)/(gamma_m1*m_kg))
        M = np.array ([[np.cos (wc*dt_s), -np.sin (wc*dt_s)],[np.sin (wc*dt_s),np.cos (wc*dt_s)]])
        v2[1],v2[2] = M@v1[[1,2]]
        
        
        gamma_m2 =  vel2gamma(v2)
        dv_dt_m0s2 = (q_C/(gamma_m2*m_kg))*  (E_V0m  - v2*np.dot(E_V0m, v2) /c_m0s**2)
        v_new_m0s = v2 + dv_dt_m0s2*dt_s/2
        r_new_mm = (R_m_i + v_new_m0s*dt_s)*1e3
        gamma = vel2gamma(v_new_m0s)
        

        v_minus_half_vec_m0s.append(v_new_m0s)
        R_vec_mm.append(r_new_mm)
        gamma_vec.append(gamma)
        
    return R_vec_mm,v_minus_half_vec_m0s,gamma_vec
        


