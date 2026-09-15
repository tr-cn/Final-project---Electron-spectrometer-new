import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from types import SimpleNamespace

#plt.rcParams.update({
#   'font.family': 'serif',
#  'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
# 'mathtext.fontset': 'stix'  # Matches Times New Roman math styles seamlessly
#})

from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from spectrometer.physics import get_lorentz_acceleration, get_electric_acceleration, get_magnetic_rotation
from spectrometer.analitic_fields import *#Magenetic_field_Analitic, Electric_field_Analitic
from spectrometer.numeric_fields import *

def MeV2m0s(q_eng_MeV, m_kg):
    q_eng_J = abs(q_eng_MeV)*1e6 * 1.602*1e-19
    c_m0s = 299792458
    v0_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
    return v0_m0s

def vel2gamma (v_m0s):
    c_m0s = 299792458
    v0_mag_m0s = np.linalg.norm(v_m0s)
    return 1 / np.sqrt ( 1 - (v0_mag_m0s/c_m0s)**2 )




# def get_magnetic_field (B_type,B_T,R_current): # Need to think how to do it correctly
#     # B_type can get: Constant, analitic, numeric
#     if B_type == "const": # here B_T is value
#         return B_T
#     if B_type == "analitic": # here B_T is the class: Magenetic_field_Analitic
#         Bx, By, Bz = B_T._get_magnetic_field(R_current)
#         return Bx, By, Bz
#     if B_type == "numeric": # Here B_T is a grid of valuse that need to be interpulated for the exact positon
#         return None


# def get_electric_field (E_type, dx_mm, dy_mm, dz_mm): # Need to think how to do it correctly
   
#     return






class Integrators():
    def __init__(self,params):
        # self.experiment = experiment
        # self.q_eng_MeV = self.experiment.q_eng_MeV
        # self.m_kg = self.experiment.m_kg
        # self.q_C = self.experiment.q_C
        # self.height_mm = self.experiment.height_mm
        # self.Bx0_T = self.experiment.Bx0_T
        # self.q_eng_J = self.experiment.q_eng_MeV*1e6 * 1.602*1e-19
        # self.h_m = self.experiment.height_mm*1e-3 /2
        # self.w_m = self.experiment.width_mm*1e-3  /2
        # self.d_m = self.experiment.depth_mm*1e-3  /2
        self.__dict__.update(params)
        self.k = 0
        if self.solution == "Analitic field":
            self.mag = Magenetic_field_Analitic(Bx0_T=self.Bx0_T, width_mm=self.width_mm, depth_mm=self.depth_mm,
                                           k=self.k, fringe = self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)
            
            
            if any (self.Exz0_Vm):
                if self.Exz0_Vm[0] != 0:
                      self.Ele = Electric_field_Analitic(E0_Vm = self.Exz0_Vm[0], width_mm=self.width_mm, depth_mm=self.depth_mm, 
                                                         k=self.k, fringe=self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)  
                elif self.Exz0_Vm[2] != 0:
                    self.Ele = Electric_field_Analitic(E0_Vm = self.Exz0_Vm[2], width_mm=self.width_mm, depth_mm=self.depth_mm, 
                                                       k=self.k, fringe=self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)
                elif all(self.Exz0_Vm==0):
                    self.Ele = self.Exz0_Vm
            
            
            
        elif self.solution == "Numeric Field":
            self.psi = MagneticPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, B0_T=self.Bx0_T,
                                    pole_y_start_mm=0, pole_y_end_mm=self.depth_mm)
            self.psi._solve_potential()
            
            self.mag = MagneticField(self.psi)
            self.mag._solve_field()
            
            
            if any (self.Exz0_Vm):
                if self.Exz0_Vm[0] != 0:
                    self.phi = ElectricPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, E0_Vm=self.Exz0_Vm[0],
                                       electrode_y_start_mm=0, electrode_y_end_mm=self.depth_mm) 
                    
                elif self.Exz0_Vm[2] != 0:
                    self.phi = ElectricPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_Z_mm=self.exp_range_X_mm, E0_Vm=self.Exz0_Vm[2],
                                       electrode_y_start_mm=0, electrode_y_end_mm=self.depth_mm) 
                
                self.phi._solve_potential()
                self.Ele = ElectricField(self.phi)
                self.Ele._solve_field()


    def _analitic_sol_vel2dist (self):
        
        # An explanation of how radius and the velocities are calculated is given
        # in the documation

        R_m = self.gamma0 * self.m_kg*self.v0_m0s / (self.q_C * self.Bx0_T)
        # R_mm = R_m *1e3
        signs = np.sign(R_m) 
        
        R_m = abs(R_m)
        Y_mm = np.sqrt ( 2*R_m[1]*self.h_m - self.h_m**2 )*1e3 #* signs[1]
        
        return Y_mm


    
    def _is_in_spectrometer(self,R_current_m):#, h_m, d_m, w_m, shield_mm, pinhole_dia_mm):
        
        if  R_current_m[1]<=0:
            
            return abs(R_current_m[2]) < self.pinhole_rad_m  and abs(R_current_m[0]) < self.pinhole_rad_m
        
        return abs(R_current_m[2]) < self.h_m and R_current_m[1] < self.d_m and abs(R_current_m[0]) < abs (self.w_m)
    
    
    def _get_magnetic_field(self,R_current_m):
        if self.solution == "Analitic field":
            if self.fringe == 0:
                if R_current_m[1]<0:
                    return np.array([0,0,0])
                else:
                    return np.array([self.Bx0_T,0,0])
                
            if self.fringe == 1:
                # mag = Magenetic_field_Analitic(self.Bx0_T, self.width_mm, self.depth_mm, self.k, self.fringe, self.sharp_edge, self.yoke)
                return  self.mag._get_magnetic_field(R_current_m)
        elif self.solution ==  "Numeric Field":
               # print(R_current_m*1e3)
               B_p =  self.mag._get_magnetic_field(R_current_m)
               
               return B_p
            
        
            
    def _get_electric_field(self,R_current_m):
        if self.solution == "Analitic field":
            if self.fringe == 0:
                if R_current_m[1]<0:
                    return [0,0,0]
                else:
                    return self.Exz0_Vm
                
            if self.fringe == 1:
                
                if any (self.Exz0_Vm):
                    if self.Exz0_Vm[0] != 0:
                        E_p = self.Ele._get_electric_field(R_current_m)
                        return E_p
                    elif self.Exz0_Vm[2] != 0:
                        E_p = self.Ele._get_electric_field(R_current_m)
                        E_pp = np.array([E_p[2],E_p[1],E_p[0]])
                        return E_pp
                else:
                    return np.array([0,0,0])
                        
                # mag = Magenetic_field_Analitic(self.Bx0_T, self.width_mm, self.depth_mm, self.k, self.fringe, self.sharp_edge, self.yoke)
                # return  mag._get_magnetic_field(R_current_m)
                
                return None
                
        elif self.solution ==  "Numeric Field":
               # print(R_current_m*1e3)
               E_p =  self.Ele._get_electric_field(R_current_m)
               
               return E_p

    def _euler (self):
        
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = self.v0_m0s#np.copy(self.v0_m0s)
        
        self.v_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);
        
        #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
        # T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
        # dt_s = T_cyclotron/steps
        
        dt_s = self.dt_s
        while self._is_in_spectrometer(R_current_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
            B_T = self._get_magnetic_field(R_current_m);
            # print(R_current_m)
            # print(B_T)
            # print('')
            E_Vm = self._get_electric_field(R_current_m)
            
            gamma = gamma_current
            v_m0s = v_current_m0s
            R_m = R_current_m        
            dv_dt_m0s2= get_lorentz_acceleration(self.q_C,gamma,self.m_kg,E_Vm,v_m0s,B_T)
            
            v_current_m0s = v_m0s + dv_dt_m0s2*dt_s
            
            dr_m = v_current_m0s*dt_s;
            R_current_m = R_m + dr_m 
            gamma_current = vel2gamma(v_current_m0s)
            
            
            self.v_vec_m0s.append(v_current_m0s)
            self.R_vec_m.append(R_current_m)
            self.gamma_vec.append(gamma_current)
            
        
        frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        gamma_exit = vel2gamma(v_exit_m0s)
        
        self.R_vec_m[-1] = R_exit_m
        self.v_vec_m0s[-1] = v_exit_m0s
        self.gamma_vec[-1] = gamma_exit
        
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec
        

    def _RK2 (self):
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = self.v0_m0s#np.copy(self.v0_m0s)
        
        self.v_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);
        
        #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
        # T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
        # dt_s = T_cyclotron/steps
        
        dt_s = self.dt_s
        while self._is_in_spectrometer(R_current_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
            B_T = self._get_magnetic_field(R_current_m);
            # print(R_current_m)
            # print(B_T)
            # print('')
            E_Vm = self._get_electric_field(R_current_m)
            
          
            v_m0s_i = v_current_m0s; 
            gamma_i = gamma_current
            R_m_i =  R_current_m 
            
            dv_dt_m0s2_i= get_lorentz_acceleration(self.q_C,gamma_i,self.m_kg,E_Vm,v_m0s_i,B_T)
            dr_dt_m0s_i = v_m0s_i;
            
            k1_v = 1/2*dv_dt_m0s2_i*dt_s
            k1_r = 1/2*dr_dt_m0s_i*dt_s
            
            
            v_m0s_m = v_m0s_i + k1_v; 
            gamma_m = vel2gamma(v_m0s_m)
            R_m_m = R_m_i + k1_r; # not realy neaded
            
            B_T_m =B_T = self._get_magnetic_field(R_m_m);
            E_Vm_m = self._get_electric_field(R_m_m)
            
            dv_dt_m0s2_m = get_lorentz_acceleration(self.q_C,gamma_m,self.m_kg,E_Vm_m,v_m0s_m,B_T_m)
            dr_dt_m0s_m = v_m0s_m;
            
            
            k2_v = dt_s * dv_dt_m0s2_m
            k2_r = dt_s * dr_dt_m0s_m
            
            
            v_current_m0s = v_m0s_i + k2_v
            R_current_m = (R_m_i + k2_r)
            gamma_current = vel2gamma(v_current_m0s)
            
            
            self.v_vec_m0s.append(v_current_m0s)
            self.R_vec_m.append(R_current_m)
            self.gamma_vec.append(gamma_current)
            
        frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        gamma_exit = vel2gamma(v_exit_m0s)
        
        self.R_vec_m[-1] = R_exit_m
        self.v_vec_m0s[-1] = v_exit_m0s
        self.gamma_vec[-1] = gamma_exit
        
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec
        

    def _RK4 (self):
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = self.v0_m0s#np.copy(self.v0_m0s)
        
        self.v_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);
        
        #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
        # T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
        # dt_s = T_cyclotron/steps
        
        dt_s = self.dt_s
        while self._is_in_spectrometer(R_current_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
            B_T = self._get_magnetic_field(R_current_m);
            # print(R_current_m)
            # print(B_T)
            # print('')
            E_Vm = self._get_electric_field(R_current_m)
            
          
            v_m0s_i = v_current_m0s
            gamma_i = gamma_current
            R_m_i =  R_current_m 
            
            dv_dt_m0s2_i = get_lorentz_acceleration(self.q_C,gamma_i,self.m_kg,E_Vm,v_m0s_i,B_T)
            dr_dt_m0s_i = v_m0s_i;
            
            k1_v = dv_dt_m0s2_i*dt_s
            k1_r = dr_dt_m0s_i*dt_s
            
                    
            v_m0s_m1 = v_m0s_i + k1_v/2; 
            gamma_m1= vel2gamma(v_m0s_m1)
            R_m_m1 = R_m_i + k1_r/2; # not realy neaded
            
            B_T_m1 = self._get_magnetic_field(R_m_m1)
            E_Vm_m1 = self._get_electric_field(R_m_m1)
            
            dv_dt_m0s2_m1 = get_lorentz_acceleration(self.q_C,gamma_m1,self.m_kg,E_Vm_m1,v_m0s_m1,B_T_m1)
            dr_dt_m0s_m1 = v_m0s_m1;        
                 
            
            k2_v = dv_dt_m0s2_m1*dt_s
            k2_r = dr_dt_m0s_m1*dt_s
            
            v_m0s_m2 = v_m0s_i + k2_v/2; 
            gamma_m2= vel2gamma(v_m0s_m2)
            R_m_m2 = R_m_i + k2_r/2; # not realy neaded
            
            B_T_m2 = self._get_magnetic_field(R_m_m2);
            E_Vm_m2 = self._get_electric_field(R_m_m2)
            
            
            dv_dt_m0s2_m2 = get_lorentz_acceleration(self.q_C,gamma_m2,self.m_kg,E_Vm_m2,v_m0s_m2,B_T_m2)
            dr_dt_m0s_m2 = v_m0s_m2;   
            
            
            k3_v = dv_dt_m0s2_m2*dt_s
            k3_r = dr_dt_m0s_m2*dt_s
            
            
            v_m0s_m3 = v_m0s_i + k3_v; 
            gamma_m3= vel2gamma(v_m0s_m3)
            R_m_m3 = R_m_i + k3_r; # not realy neaded
            
            B_T_m3 = self._get_magnetic_field(R_m_m3)
            E_Vm_m3 = self._get_electric_field(R_m_m3)
            
            dv_dt_m0s2_m3 = get_lorentz_acceleration(self.q_C,gamma_m3,self.m_kg,E_Vm_m3,v_m0s_m3,B_T_m3)
            dr_dt_m0s_m3 = v_m0s_m3;  
            
            k4_v = dv_dt_m0s2_m3*dt_s
            k4_r = dr_dt_m0s_m3*dt_s 
            
            
            v_current_m0s = v_m0s_i + 1/6 * (k1_v + 2*k2_v + 2*k3_v + k4_v)
            R_current_m = (R_m_i + 1/6 * (k1_r + 2*k2_r + 2*k3_r + k4_r))
            gamma_current = vel2gamma(v_current_m0s)
            
            
            self.v_vec_m0s.append(v_current_m0s)
            self.R_vec_m.append(R_current_m)
            self.gamma_vec.append(gamma_current)
            
        frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        gamma_exit = vel2gamma(v_exit_m0s)
        
        self.R_vec_m[-1] = R_exit_m
        self.v_vec_m0s[-1] = v_exit_m0s
        self.gamma_vec[-1] = gamma_exit
        
        
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec
            

    def _Boris_pusher (self):
        dt_s = self.dt_s
        
        R_current_m = self.R0_m
        B_T = self._get_magnetic_field(R_current_m);
        E_Vm = self._get_electric_field(R_current_m)
        dv_dt_0_m0s2 = get_lorentz_acceleration(self.q_C,self.gamma0,self.m_kg,E_Vm,self.v0_m0s,B_T)
        
        v_minus_half_m0s = self.v0_m0s - 1/2 * dv_dt_0_m0s2 * dt_s
        
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = v_minus_half_m0s
        
        self.v_minus_half_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_minus_half_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);

    
        while self._is_in_spectrometer(R_current_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            B_T = self._get_magnetic_field(R_current_m);
            # print(R_current_m)
            # print(B_T)
            # print('')
            E_Vm_i = self._get_electric_field(R_current_m)# change it
            
            v_minus_half_m0s = v_current_m0s
            gamma_i = gamma_current
            R_m_i = R_current_m 
            
            dv1_m0s= get_electric_acceleration(self.q_C,gamma_i,self.m_kg,E_Vm_i,v_minus_half_m0s,B_T=None) * dt_s/2
            v1_m0s = v_minus_half_m0s + dv1_m0s 
            gamma_m1 =  vel2gamma(v1_m0s)
            # R_m_m1 = R_m_i + dv1_m0s*dt_s/2

            
            
           
            
            v2_m0s  = get_magnetic_rotation(self.q_C,gamma_m1,self.m_kg,E_Vm_i,v1_m0s,B_T,dt_s)
            gamma_m2 =  vel2gamma(v2_m0s)
            
            
            v_current_m0s = v2_m0s + get_electric_acceleration(self.q_C,gamma_m2,self.m_kg,E_Vm_i,v2_m0s,B_T=None)*dt_s/2
            R_current_m = (R_m_i + v_current_m0s*dt_s)
            gamma_current = vel2gamma(v_current_m0s)
            
    
            self.v_minus_half_vec_m0s.append(v_current_m0s)
            self.R_vec_m.append(R_current_m)
            
            
            self.gamma_vec.append(gamma_current)
            
        frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        gamma_exit = vel2gamma(v_exit_m0s)
        
        self.R_vec_m[-1] = R_exit_m
        self.v_vec_m0s[-1] = v_exit_m0s
        self.gamma_vec[-1] = gamma_exit
        
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_minus_half_vec_m0s,self.gamma_vec
        



