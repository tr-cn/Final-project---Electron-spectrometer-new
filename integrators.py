import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#plt.rcParams.update({
#   'font.family': 'serif',
#  'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
# 'mathtext.fontset': 'stix'  # Matches Times New Roman math styles seamlessly
#})

from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from spectrometer.physics import get_lorentz_acceleration, get_electric_acceleration, get_magnetic_rotation




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


def is_in_spectrometer(R_current_m,h_m,d_m,w_m,shield_mm,pinhole_dia_mm):
    if  R_current_m[1]<=0:
        pinhole_rad_m = pinhole_dia_mm*1e-3/2
        return abs(R_current_m[2]) < pinhole_rad_m  and abs(R_current_m[0]) < pinhole_rad_m
    
    return abs(R_current_m[2]) < h_m and R_current_m[1] < d_m and abs(R_current_m[0]) < abs (w_m)
    

def euler (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm,pinhole_dia_mm, fringe):
    h_m = height_mm/2*1e-3; d_m = depth_mm*1e-3; w_m = width_mm/2*1e-3;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg); v_current_m0s=v0_m0s;
    R0_m = R0_mm *1e-3; R_current_m= R0_m;
    
    gamma = vel2gamma(v0_m0s); gamma_current = gamma
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_m = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_m.append(R0_m);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/steps
    B_in_T = B_T;
    B_zero_T = np.array([0,0,0])
    
    while is_in_spectrometer(R_current_m,h_m,d_m,w_m,shield_mm,pinhole_dia_mm):
        
        B_T = B_in_T
        if fringe==0 and R_current_m[1]<0:
            B_T = B_zero_T  
        
        gamma = gamma_current
        v_m0s = v_current_m0s
        R_m = R_current_m        
        dv_dt_m0s2= get_lorentz_acceleration(q_C,gamma,m_kg,E_V0m,v_m0s,B_T)
        
        v_current_m0s = v_m0s + dv_dt_m0s2*dt_s
        
        dr_m = v_current_m0s*dt_s;
        R_current_m = R_m + dr_m 
        gamma_current = vel2gamma(v_current_m0s)
        
        
        v_vec_m0s.append(v_current_m0s)
        R_vec_m.append(R_current_m)
        gamma_vec.append(gamma_current)
    
    
    R_vec_mm = [R_m*1e3 for R_m in R_vec_m]
    return R_vec_mm,v_vec_m0s,gamma_vec
        

def RK2 (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe):
    h_m = height_mm/2*1e-3; d_m = depth_mm*1e-3; w_m = width_mm/2*1e-3;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg); v_current_m0s=v0_m0s;
    R0_m = R0_mm *1e-3; R_current_m= R0_m;
    
    gamma = vel2gamma(v0_m0s); gamma_current = gamma
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_m = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_m.append(R0_m);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/steps
    B_in_T = B_T;
    B_zero_T = np.array([0,0,0])
    while is_in_spectrometer(R_current_m,h_m,d_m,w_m,shield_mm,pinhole_dia_mm):
        
        B_T = B_in_T
        if fringe==0 and R_current_m[1]<0:
            B_T = B_zero_T 
        
      
        v_m0s_i = v_current_m0s; 
        gamma_i = gamma_current
        R_m_i =  R_current_m 
        
        dv_dt_m0s2_i= get_lorentz_acceleration(q_C,gamma_i,m_kg,E_V0m,v_m0s_i,B_T)
        dr_dt_m0s_i = v_m0s_i;
        
        k1_v = 1/2*dv_dt_m0s2_i*dt_s
        k1_r = 1/2*dr_dt_m0s_i*dt_s
        
        
        v_m0s_m = v_m0s_i + k1_v; 
        gamma_m = vel2gamma(v_m0s_m)
        R_m_m = R_m_i + k1_r; # not realy neaded
        
        
        dv_dt_m0s2_m = get_lorentz_acceleration(q_C,gamma_m,m_kg,E_V0m,v_m0s_m,B_T)
        dr_dt_m0s_m = v_m0s_m;
        
        
        k2_v = dt_s * dv_dt_m0s2_m
        k2_r = dt_s * dr_dt_m0s_m
        
        
        v_current_m0s = v_m0s_i + k2_v
        R_current_m = (R_m_i + k2_r)
        gamma_current = vel2gamma(v_current_m0s)
        
        
        v_vec_m0s.append(v_current_m0s)
        R_vec_m.append(R_current_m)
        gamma_vec.append(gamma_current)
        
    R_vec_mm = [R_m*1e3 for R_m in R_vec_m]
    return R_vec_mm,v_vec_m0s,gamma_vec
        

def RK4 (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm,pinhole_dia_mm, fringe):
    h_m = height_mm/2*1e-3; d_m = depth_mm*1e-3; w_m = width_mm/2*1e-3;
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg); v_current_m0s=v0_m0s;
    R0_m = R0_mm *1e-3; R_current_m= R0_m;
    
    gamma = vel2gamma(v0_m0s); gamma_current = gamma
    
    v_vec_m0s = list([]); gamma_vec = list([]); R_vec_m = list([]); 
    v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_m.append(R0_m);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
    dt_s = T_cyclotron/steps
    B_in_T = B_T;
    B_zero_T = np.array([0,0,0])
    while is_in_spectrometer(R_current_m,h_m,d_m,w_m,shield_mm,pinhole_dia_mm):
        
        B_T = B_in_T
        if fringe==0 and R_current_m[1]<0:
            B_T = B_zero_T 
        
      
        v_m0s_i = v_current_m0s; 
        gamma_i = gamma_current
        R_m_i =  R_current_m 
        
        dv_dt_m0s2_i = get_lorentz_acceleration(q_C,gamma_i,m_kg,E_V0m,v_m0s_i,B_T)
        dr_dt_m0s_i = v_m0s_i;
        
        k1_v = dv_dt_m0s2_i*dt_s
        k1_r = dr_dt_m0s_i*dt_s
        
                
        v_m0s_m1 = v_m0s_i + k1_v/2; 
        gamma_m1= vel2gamma(v_m0s_m1)
        R_m_m1 = R_m_i + k1_r/2; # not realy neaded
        
        dv_dt_m0s2_m1 = get_lorentz_acceleration(q_C,gamma_m1,m_kg,E_V0m,v_m0s_m1,B_T)
        dr_dt_m0s_m1 = v_m0s_m1;        
             
        
        k2_v = dv_dt_m0s2_m1*dt_s
        k2_r = dr_dt_m0s_m1*dt_s
        
        v_m0s_m2 = v_m0s_i + k2_v/2; 
        gamma_m2= vel2gamma(v_m0s_m2)
        R_m_m2 = R_m_i + k2_r/2; # not realy neaded
        
        
        dv_dt_m0s2_m2 = get_lorentz_acceleration(q_C,gamma_m2,m_kg,E_V0m,v_m0s_m2,B_T)
        dr_dt_m0s_m2 = v_m0s_m2;   
        
        
        k3_v = dv_dt_m0s2_m2*dt_s
        k3_r = dr_dt_m0s_m2*dt_s
        
        
        v_m0s_m3 = v_m0s_i + k3_v; 
        gamma_m3= vel2gamma(v_m0s_m3)
        R_m_m3 = R_m_i + k3_r; # not realy neaded
        
        dv_dt_m0s2_m3 = get_lorentz_acceleration(q_C,gamma_m3,m_kg,E_V0m,v_m0s_m3,B_T)
        dr_dt_m0s_m3 = v_m0s_m3;  
        
        k4_v = dv_dt_m0s2_m3*dt_s
        k4_r = dr_dt_m0s_m3*dt_s 
        
        
        v_current_m0s = v_m0s_i + 1/6 * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        R_current_m = (R_m_i + 1/6 * (k1_r + 2*k2_r + 2*k3_r + k4_r))
        gamma_current = vel2gamma(v_current_m0s)
        
        
        v_vec_m0s.append(v_current_m0s)
        R_vec_m.append(R_current_m)
        gamma_vec.append(gamma_current)
    
    R_vec_mm = [R_m*1e3 for R_m in R_vec_m]
    
    return R_vec_mm,v_vec_m0s,gamma_vec
        

def Boris_pusher (q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm,pinhole_dia_mm, fringe):
    
    h_m = height_mm/2*1e-3; d_m = depth_mm*1e-3; w_m = width_mm/2*1e-3;
    
    gamma_0 = vel2gamma( MeV2m0s(q_eng_MeV, m_kg))
    T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma_0*m_kg) )**-1
    dt_s = T_cyclotron/steps
    
    
    B_in_T = B_T;
    B_zero_T = np.array([0,0,0])
    
    if fringe==0 and R0_mm[1]*1e-3<0:
        B_T = B_zero_T
    
    v0_m0s = MeV2m0s(q_eng_MeV, m_kg)
    dv_dt_0_m0s2 = get_lorentz_acceleration(q_C,gamma_0,m_kg,E_V0m,v0_m0s,B_T)
    v_minus_half_m0s = v0_m0s - 1/2 * dv_dt_0_m0s2 * dt_s
    
    
    v_current_m0s=v_minus_half_m0s
    R0_m = R0_mm *1e-3; R_current_m= R0_m;
    gamma = vel2gamma(v_current_m0s); gamma_current = gamma
    

    
    
    v_minus_half_vec_m0s = list([]); gamma_vec = list([]); R_vec_m = list([]); 
    v_minus_half_vec_m0s.append(v_minus_half_m0s); gamma_vec.append(gamma); R_vec_m.append(R0_m);
    
    #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
    
    
       
    

    while is_in_spectrometer(R_current_m,h_m,d_m,w_m,shield_mm,pinhole_dia_mm):
        
        B_T = B_in_T
        if fringe==0 and R_current_m[1]<0:
            B_T = B_zero_T
        v_minus_half_m0s = v_current_m0s
        gamma_i = gamma_current
        R_m_i = R_current_m 
        
        v1_m0s = v_minus_half_m0s + get_electric_acceleration(q_C,gamma_i,m_kg,E_V0m,v_minus_half_m0s,B_T) * dt_s/2
        gamma_m1 =  vel2gamma(v1_m0s)
        
        
        v2_m0s  = get_magnetic_rotation(q_C,gamma_m1,m_kg,E_V0m,v1_m0s,B_T,dt_s)
        gamma_m2 =  vel2gamma(v2_m0s)
        
        v_current_m0s = v2_m0s + get_electric_acceleration(q_C,gamma_m2,m_kg,E_V0m,v2_m0s,B_T)*dt_s/2
        R_current_m = (R_m_i + v_current_m0s*dt_s)
        gamma_current = vel2gamma(v_current_m0s)
        

        v_minus_half_vec_m0s.append(v_current_m0s)
        R_vec_m.append(R_current_m)
        
        
        gamma_vec.append(gamma_current)
    R_vec_mm = [R_m*1e3 for R_m in R_vec_m]
    return R_vec_mm,v_minus_half_vec_m0s,gamma_vec
        



