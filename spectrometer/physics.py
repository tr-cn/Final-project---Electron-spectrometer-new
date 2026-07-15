import numpy as np



def get_lorentz_acceleration(q_C,gamma,m_kg,E_V0m,v_m0s,B_T):
    c_m0s = 299792458
    
    dv_dt_m0s2= (q_C/(gamma*m_kg)) * (E_V0m + np.cross(v_m0s,B_T) -  v_m0s*np.dot(E_V0m, v_m0s) /c_m0s**2 )
    
    return dv_dt_m0s2


def get_electric_acceleration(q_C,gamma,m_kg,E_V0m,v_m0s,B_T):
    c_m0s = 299792458
    
    return (q_C/(gamma*m_kg)) * (E_V0m  -  v_m0s*np.dot(E_V0m, v_m0s) /c_m0s**2 )
    

def get_magnetic_rotation(q_C,gamma_m1,m_kg,E_V0m,v1_m0s,B_T,dt_s):
    c_m0s = 299792458
    
    t_vec = (q_C / (m_kg * gamma_m1)) * B_T * (dt_s / 2)  # for general boris algorithm 
    t_mag2 = np.dot(t_vec, t_vec)
    s_vec = 2 * t_vec / (1 + t_mag2)                        
        
    v_prime_m0s = v1_m0s + np.cross(v1_m0s,t_vec)
    v2_m0s  = v1_m0s + np.cross(v_prime_m0s,s_vec)
    return v2_m0s