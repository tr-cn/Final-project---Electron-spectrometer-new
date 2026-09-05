import matplotlib.pyplot as plt
import numpy as np
import time


from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from spectrometer.physics import get_lorentz_acceleration, get_electric_acceleration, get_magnetic_rotation
from spectrometer.fields import Magenetic_field_Analitic




from spectrometer.geometry import geometry
import spectrometer.integrators as integ
import spectrometer.plot_result as pr
import spectrometer.fields as fields




def Magnetic_Field_Init(B_type, B0_T):


def get_magnetic_field (B_type,B_T,R_current): # Need to think how to do it correctly
    # B_type can get: Constant, analitic, numeric
    if B_type == "const": # here B_T is value
        return B_T
    if B_type == "analitic": # here B_T is the class: Magenetic_field_Analitic
        Bx, By, Bz = B_T._get_magnetic_field(R_current)
        return Bx, By, Bz
    if B_type == "numeric": # Here B_T is a grid of valuse that need to be interpulated for the exact positon
        return None





class Integrators():
    def __init__(self, q_eng_MeV, m_kg, q_C, height_mm, width_mm, depth_mm, B_type, B0_T, E_V0m, R0_mm, steps, shield_mm,pinhole_dia_mm, fringe):
        self.q_eng_MeV = q_eng_MeV
        self.m_kg = m_kg 
        self.q_C = q_C
        self.height_m = height_mm*1e-3
        self.h_m = self.height_m/2
        self.width_m =  width_mm*1e-3
        self.w_m = self.width_m/2
        self.d_m = depth_mm*1e-3
        self.B_type = B_type
        self.B0_T = B0_T
        self.E_V0m = E_V0m
        self.R0_m = R0_mm*1e-3
        self.steps = steps
        self.shield_m = shield_mm*1e-3
        self.pinhole_rad_m = 0.5 * pinhole_dia_mm*1e-3 
        self.fring = fringe
        self.R_current_m = None
        
        
        def _MeV2m0s(self):
            q_eng_MeV = self.q_eng_MeV; m_kg= self.m_kg;
            q_eng_J = abs(q_eng_MeV)*1e6 * 1.602*1e-19
            c_m0s = 299792458
            v0_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
            return v0_m0s

        def _vel2gamma (self,v_m0s):
            c_m0s = 299792458
            v0_mag_m0s = np.linalg.norm(v_m0s)
            return 1 / np.sqrt ( 1 - (v0_mag_m0s/c_m0s)**2 )
        
        
        
        def _is_in_spectrometer(self,R_current_m):
            R_current_m
            if  self.R_current_m[1]<=0:
                return abs(R_current_m[2]) < self.pinhole_rad_m  and abs(self.R_current_m[0]) < self.pinhole_rad_m
            
            return abs(R_current_m[2]) < self.h_m and R_current_m[1] < self.d_m and abs(R_current_m[0]) < abs (self.w_m)
        
        
        def __euler (self):
                     
            v0_m0s = self._MeV2m0s(); 
            v_current_m0s = v0_m0s;
            R_current_m = self.R0_m;
            
            gamma = self._vel2gamma(v0_m0s); gamma_current = gamma
            
            v_vec_m0s = list([]); gamma_vec = list([]); R_vec_m = list([]); 
            v_vec_m0s.append(v0_m0s); gamma_vec.append(gamma); R_vec_m.append(self.R0_m);
            
            #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
            T_cyclotron = ( abs(q_C)*np.linalg.norm(self.B0_T)/(np.pi*gamma*m_kg) )**-1
            dt_s = T_cyclotron/steps
            
            B_in_T = get_magnetic_field()            
            B_in_T = self.B0_T;
            B_zero_T = np.array([0,0,0])
            
            while _is_in_spectrometer(R_current_m):
                
                B0_T = B_in_T
                if fringe==0 and R_current_m[1]<0:
                    B0_T = B_zero_T  
                
                gamma = gamma_current
                v_m0s = v_current_m0s
                R_m = R_current_m        
                dv_dt_m0s2= get_lorentz_acceleration(q_C,gamma,m_kg,E_V0m,v_m0s,B0_T)
                
                v_current_m0s = v_m0s + dv_dt_m0s2*dt_s
                
                dr_m = v_current_m0s*dt_s;
                R_current_m = R_m + dr_m 
                gamma_current = _vel2gamma(v_current_m0s)
                
                
                v_vec_m0s.append(v_current_m0s)
                R_vec_m.append(R_current_m)
                gamma_vec.append(gamma_current)
            
            
            R_vec_mm = [R_m*1e3 for R_m in R_vec_m]
            return R_vec_mm,v_vec_m0s,gamma_vec





if __name__ == "__main__":
    plt.close('all')
    me_kg = 9.109*1e-31
    B_type = "ana"
    B0_T = np.array([0.5,0.0,0.0])
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

