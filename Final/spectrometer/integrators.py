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
        
            if self.Ex0_Vm != 0:
                self.Ele = Electric_field_Analitic(Ex0_Vm = self.Ex0_Vm, width_mm=self.width_mm, depth_mm=self.depth_mm, 
                                                             k=self.k, fringe=self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)              
               
        
        
        elif self.solution == "Numeric Field":
            self.psi = MagneticPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, B0_T=self.Bx0_T,
                                    pole_y_start_mm=0, pole_y_end_mm=self.depth_mm)
            self.psi._solve_potential()
            
            self.mag = MagneticField(self.psi, self.grid_interpulation)
            self.mag._solve_field()
            
            
            if self.Ex0_Vm != 0:
                self.phi = ElectricPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, Ex0_Vm=self.Ex0_Vm,
                                       electrode_y_start_mm=0, electrode_y_end_mm=self.depth_mm) 
                self.phi._solve_potential()
                self.Ele = ElectricField(self.phi, self.grid_interpulation)
                self.Ele._solve_field()
                

        
    def _update_params(self,params):
        self.__dict__.update(params)
        self.params = params
        
        
    
    def _update_solution(self,solution):
        self.solution = solution
        if self.solution == "Analitic field":
            self.mag = Magenetic_field_Analitic(Bx0_T=self.Bx0_T, width_mm=self.width_mm, depth_mm=self.depth_mm,
                                           k=self.k, fringe = self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)
        
            if self.Ex0_Vm != 0:
                self.Ele = Electric_field_Analitic(Ex0_Vm = self.Ex0_Vm, width_mm=self.width_mm, depth_mm=self.depth_mm, 
                                                             k=self.k, fringe=self.fringe, sharp_edge=self.sharp_edge, yoke=self.yoke)              
               
        
        
        elif self.solution == "Numeric Field":
            self.psi = MagneticPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, B0_T=self.Bx0_T,
                                    pole_y_start_mm=0, pole_y_end_mm=self.depth_mm)
            self.psi._solve_potential()
            
            self.mag = MagneticField(self.psi,self.grid_interpulation)
            self.mag._solve_field()
            
            
            if self.Ex0_Vm != 0:
                self.phi = ElectricPotential(Ny=self.Ny_p, Nx=self.Nx_p, exp_range_Y_mm=self.exp_range_Y_mm, exp_range_X_mm=self.exp_range_X_mm, Ex0_Vm=self.Ex0_Vm,
                                       electrode_y_start_mm=0, electrode_y_end_mm=self.depth_mm) 
                self.phi._solve_potential()
                self.Ele = ElectricField(self.phi,self.grid_interpulation)
                self.Ele._solve_field()

    def _analitic_sol_vel2dist (self):
        
        # An explanation of how radius and the velocities are calculated is given
        # in the documation

        R_m = self.gamma0 * self.m_kg*self.v0_m0s / (self.q_C * self.Bx0_T)
        # R_mm = R_m *1e3
        signs = np.sign(R_m) 
        
        R_m = abs(R_m[1])
        h = self.h_m - self.R0_m[2]
        Y_mm = np.sqrt ( 2*R_m*h - h**2 )*1e3 
        theta = np.arccos( (R_m-h)/ R_m)
        tau = theta*self.m_kg/(self.q_C*self.Bx0_T)
        
        c_m0s = 299792458
        
        kappa = self.Ex0_Vm*self.q_C/(self.m_kg * c_m0s)
        X_m = self.R0_m[0] +  self.gamma0*self.v0_m0s[0]*np.sinh(kappa*tau) + c_m0s*self.gamma0/kappa*(np.cosh(kappa*tau)-1 )
        X_mm = X_m*1e3
        
        return X_mm, Y_mm


    
    def _is_in_spectrometer(self,R_current_m,R_prev_m):#, h_m, d_m, w_m, shield_mm, pinhole_dia_mm):
        
        if  R_current_m[1]<=0:
            
            return abs(R_current_m[2]) < self.pinhole_rad_m  and abs(R_current_m[0]) < self.pinhole_rad_m and R_current_m[1]>R_prev_m[1]
        
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
               if self.fringe == 0:
                   if R_current_m[1]<0:
                       return np.array([0,0,0])
             
               B_p =  self.mag._get_magnetic_field(R_current_m)
               
               return B_p
            
        
            
    def _get_electric_field(self,R_current_m):
        if self.solution == "Analitic field":
            if self.fringe == 0:
                if R_current_m[1]<0:
                    return [0,0,0]
                else:
                    return np.array([self.Ex0_Vm,0,0])
                
            if self.fringe == 1:
                
                if self.Ex0_Vm != 0:
                    E_p = self.Ele._get_electric_field(R_current_m)
                else:
                    E_p = np.array([0,0,0])
                return E_p
                        
                # mag = Magenetic_field_Analitic(self.Bx0_T, self.width_mm, self.depth_mm, self.k, self.fringe, self.sharp_edge, self.yoke)
                # return  mag._get_magnetic_field(R_current_m)
                
                return None
                
        elif self.solution ==  "Numeric Field":
               # print(R_current_m*1e3)
               if self.Ex0_Vm != 0: 
                   E_p =  self.Ele._get_electric_field(R_current_m)
               else:
                   E_p = np.array([0,0,0])
               return E_p
           
    def _hermite_exit_crossing(self, R_m, v_ms):
        """
        Hermite cubic interpolation to the spectrometer boundary.
    
        Input:
            R_m  : list/array of particle positions [m]
            v_ms : list/array of particle velocities [m/s]
    
        Returns:
            R_exit_m    : exit position [m]
            v_exit_m0s  : exit velocity [m/s]
            gamma_exit  : Lorentz gamma at the exit
        """
    
        # Two endpoints of the final integration step
        R0 = np.asarray(R_m[-2], dtype=float)
        R1 = np.asarray(R_m[-1], dtype=float)
    
        v0 = np.asarray(v_ms[-2], dtype=float)
        v1 = np.asarray(v_ms[-1], dtype=float)
    
        dt = self.dt_s
        
        z_limit= self.pinhole_rad_m if R1[1]<= 0 else self.h_m
        x_limit= self.pinhole_rad_m if R1[1]<= 0 else self.w_m
    
        # Determine which physical boundary was crossed
        if R1[2] > z_limit:
            axis = 2
            target = z_limit
    
        elif R1[2] < -z_limit:
            axis = 2
            target = -z_limit
    
        elif R1[0] > x_limit:
            axis = 0
            target = x_limit
    
        elif R1[0] < -x_limit:
            axis = 0
            target = -x_limit
    
        elif R1[1] > self.d_m:
            axis = 1
            target = self.d_m
    
        elif R1[1] < 0.0:
            axis = 1
            target = 0.0
    
        else:
            # No boundary crossing: return the final ordinary RK4 point
            R_exit_m = R1.copy()
            v_exit_m0s = v1.copy()
            gamma_exit = vel2gamma(v_exit_m0s)
    
            return R_exit_m, v_exit_m0s, gamma_exit
    
        # Scalar Hermite polynomial of the coordinate that crosses the boundary
        p0 = R0[axis]
        p1 = R1[axis]
    
        # Tangents with respect to normalized coordinate s in [0,1]
        m0 = v0[axis] * dt
        m1 = v1[axis] * dt
    
        def h_basis(s):
            h00 = 2.0*s**3 - 3.0*s**2 + 1.0
            h10 = s**3 - 2.0*s**2 + s
            h01 = -2.0*s**3 + 3.0*s**2
            h11 = s**3 - s**2
    
            return h00, h10, h01, h11
    
        def crossing_coordinate(s):
            h00, h10, h01, h11 = h_basis(s)
    
            return h00*p0 + h10*m0 + h01*p1 + h11*m1
    
        # Solve R_axis(s) = target by bisection
        lo = 0.0
        hi = 1.0
    
        f_lo = crossing_coordinate(lo) - target
        f_hi = crossing_coordinate(hi) - target
    
        # Numerical safety: should not occur if this function is called only
        # after the particle has crossed a boundary.
        if f_lo * f_hi > 0.0:
            frac = (target - p0) / (p1 - p0)
            frac = np.clip(frac, 0.0, 1.0)
            t_frac = frac
    
        else:
            for _ in range(60):
                mid = 0.5*(lo + hi)
                f_mid = crossing_coordinate(mid) - target
    
                if f_lo * f_mid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    f_lo = f_mid
    
            t_frac = 0.5*(lo + hi)
    
        # Hermite position interpolation for the full vector R(s)
        h00, h10, h01, h11 = h_basis(t_frac)
    
        R_exit_m = (
            h00 * R0
            + h10 * (v0 * dt)
            + h01 * R1
            + h11 * (v1 * dt)
        )
    
        # Derivative of the Hermite trajectory:
        # v_H(s) = dR/dt = (1/dt) dR/ds
        s = t_frac
    
        dh00 = 6.0*s**2 - 6.0*s
        dh10 = 3.0*s**2 - 4.0*s + 1.0
        dh01 = -6.0*s**2 + 6.0*s
        dh11 = 3.0*s**2 - 2.0*s
    
        v_hermite_raw = (
            (dh00 / dt) * R0
            + dh10 * v0
            + (dh01 / dt) * R1
            + dh11 * v1
        )
    
        # Keep the speed interpolation consistent with _exect_exist.
        # In a magnetic-only run, |v0| and |v1| should be nearly equal.
        mag0 = np.linalg.norm(v0)
        mag1 = np.linalg.norm(v1)
    
        mag_exit = mag0 + t_frac*(mag1 - mag0)
    
        # Tangent direction from the Hermite curve
        tangent_norm = np.linalg.norm(v_hermite_raw)
    
        if tangent_norm == 0.0:
            # Extremely defensive fallback; normally never reached.
            v_exit_m0s = v0 + t_frac*(v1 - v0)
    
        else:
            dir_exit = v_hermite_raw / tangent_norm
            v_exit_m0s = mag_exit * dir_exit
    
        gamma_exit = vel2gamma(v_exit_m0s)
    
        return R_exit_m, v_exit_m0s, gamma_exit
           
    def _exect_exist(self,R_m,v_ms):
        frac =0 # for non entry situation
        z_limit= self.pinhole_rad_m if R_m[-1][2]<= 0 else self.h_m
        x_limit= self.pinhole_rad_m if R_m[-1][0]<= 0 else self.w_m
        
        if (R_m[-1][2])>z_limit:
            frac = (z_limit - R_m[-2][2]) / ((R_m[-1][2]) - R_m[-2][2])
        
        elif (R_m[-1][2])<-z_limit:
            frac = (-z_limit - R_m[-2][2]) / ((R_m[-1][2]) - R_m[-2][2])
        
        elif (R_m[-1][0])>x_limit:
            frac = (self.w_m - R_m[-2][0]) / (R_m[-1][0] - R_m[-2][0])  
        
        elif R_m[-1][0]<-x_limit:
            frac = (-self.w_m - R_m[-2][0]) / (R_m[-1][0] - R_m[-2][0]) 
        
        elif (R_m[-1][1])>self.d_m:
            frac = (self.d_m - R_m[-2][1]) / (R_m[-1][1] - R_m[-2][1])  
        
        elif (R_m[-1][1])<0:
            frac = (0 - R_m[-2][1]) / (R_m[-1][1] - R_m[-2][1])   
            
       
        v_minus2 = v_ms[-2]
        v_minus1 = v_ms[-1]
        
        mag_minus2 = np.linalg.norm(v_minus2)
        mag_minus1 = np.linalg.norm(v_minus1)
        mag_exit = mag_minus2 + frac * (mag_minus1 - mag_minus2)    
        
        dir_minus2 = v_minus2 / mag_minus2
        dir_minus1 = v_minus1 / mag_minus1
        dir_exit_raw = dir_minus2 + frac * (dir_minus1 - dir_minus2)
        dir_exit = dir_exit_raw / np.linalg.norm(dir_exit_raw)       
       
        R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        v_exit_m0s = dir_exit * mag_exit
        gamma_exit = vel2gamma(v_exit_m0s)
        
        # print(f"frac={frac:.4f}, |v_minus2|={mag_minus2:.4e}, |v_minus1|={mag_minus1:.4e}")
        
        
        return R_exit_m, v_exit_m0s, gamma_exit

        
    def _clipper(self,R_m)  :
        R_clamped = np.array(R_m, dtype=float).copy()
        R_clamped[0] = np.clip(R_clamped[0], -self.w_m, self.w_m)
        # R_clamped[1] = np.clip(R_clamped[1], 0, self.d_m)
        R_clamped[2] = np.clip(R_clamped[2], -self.h_m, self.h_m)
        return R_clamped
            
            
        
        
        

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
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
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
            
            R_prev_m = self.R_vec_m[-2]
            
        
        # frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        # R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        # v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        # gamma_exit = vel2gamma(v_exit_m0s)
        
        # self.R_vec_m[-1] = R_exit_m
        # self.v_vec_m0s[-1] = v_exit_m0s
        # self.gamma_vec[-1] = gamma_exit
        self.R_vec_m[-1], self.v_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_vec_m0s)
        
        # print(f"N_steps={self.N_steps}, actual_steps={len(self.R_vec_m)}")
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
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
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
            # R_m_m = self._clipper(R_m_m)
           
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
            
            R_prev_m = self.R_vec_m[-2]
            
        # frac = (self.h_m - abs(self.R_vec_m[-2][2])) / (abs(self.R_vec_m[-1][2]) - abs(self.R_vec_m[-2][2]))
        
        # R_exit_m = self.R_vec_m[-2] + frac * (self.R_vec_m[-1] - self.R_vec_m[-2])
        # v_exit_m0s =  self.v_vec_m0s[-2] + frac * ( self.v_vec_m0s[-1] -  self.v_vec_m0s[-2])
        # gamma_exit = vel2gamma(v_exit_m0s)
        
        # self.R_vec_m[-1] = R_exit_m
        # self.v_vec_m0s[-1] = v_exit_m0s
        # self.gamma_vec[-1] = gamma_exit
        self.R_vec_m[-1], self.v_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_vec_m0s)
        # self.R_vec_m[-1], self.v_vec_m0s[-1], frac =  self._hermite_exit_crossing(self.R_vec_m[-2], self.R_vec_m[-1], self.v_vec_m0s[-2], self.v_vec_m0s[-1], self.dt_s, self.h_m)
        # print(f"N_steps={self.N_steps}, actual_steps={len(self.R_vec_m)}")
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec
        

    def _RK4_Linear (self):
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = self.v0_m0s#np.copy(self.v0_m0s)
        
        self.v_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);
        
        #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
        # T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
        # dt_s = T_cyclotron/steps
        
        dt_s = self.dt_s
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
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
            
            R_prev_m = self.R_vec_m[-2]

        self.R_vec_m[-1], self.v_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_vec_m0s)
        # self.R_vec_m[-1], self.v_vec_m0s[-1], self.gamma_vec[-1] = self._hermite_exit_crossing(self.R_vec_m, self.v_vec_m0s)
        
        # print(f"N_steps={self.N_steps}, actual_steps={len(self.R_vec_m)}")
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec

    def _RK4_Hermit (self):
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = self.v0_m0s#np.copy(self.v0_m0s)
        
        self.v_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);
        
        #CFL = 0.00000001; dx_m = 1e-4; dt_s = dx_m*CFL; # not realy neaded in euler
        # T_cyclotron = ( abs(q_C)*np.linalg.norm(B_T)/(np.pi*gamma*m_kg) )**-1
        # dt_s = T_cyclotron/steps
        
        dt_s = self.dt_s
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            
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
            
            R_prev_m = self.R_vec_m[-2]

        # self.R_vec_m[-1], self.v_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_vec_m0s)
        self.R_vec_m[-1], self.v_vec_m0s[-1], self.gamma_vec[-1] = self._hermite_exit_crossing(self.R_vec_m, self.v_vec_m0s)
        # self.R_vec_m[-1], self.v_vec_m0s[-1], frac =  self._hermite_exit_crossing(self.R_vec_m[-2], self.R_vec_m[-1], self.v_vec_m0s[-2], self.v_vec_m0s[-1], self.dt_s, self.h_m)
        # print(f"N_steps={self.N_steps}, actual_steps={len(self.R_vec_m)}")
        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        
        return self.R_vec_mm,self.v_vec_m0s,self.gamma_vec
            
    def _Boris_pusher (self):
        dt_s = self.dt_s
        
        R_current_m = self.R0_m
        # B_T = self._get_magnetic_field(R_current_m);
        E_Vm = self._get_electric_field(R_current_m)
        dv_dt_0_m0s2 = get_electric_acceleration(self.q_C,self.gamma0,self.m_kg,E_Vm,self.v0_m0s,B_T=None)
        
        v_minus_half_m0s = self.v0_m0s - 1/2 * dv_dt_0_m0s2 * dt_s
        
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = v_minus_half_m0s
        
        self.v_minus_half_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_minus_half_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);

        count = 0
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            count +=1
            # if count%11681 == 0:
            #     print (count)
            #     print (R_current_m)
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
            
            R_prev_m = self.R_vec_m[-2]
            
            
        self.R_vec_m[-1], self.v_minus_half_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_minus_half_vec_m0s)
        # self.R_vec_m[-1], self.v_minus_half_vec_m0s[-1], frac =  self._hermite_exit_crossing(self.R_vec_m[-2], self.R_vec_m[-1], self.v_minus_half_vec_m0s[-2], self.v_minus_half_vec_m0s[-1], self.dt_s, self.h_m)
        

        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_minus_half_vec_m0s,self.gamma_vec
        
    
    def _Collocated_Boris_pusher (self): # https://arxiv.org/pdf/2607.12272
        dt_s = self.dt_s
        
        R_current_m = self.R0_m
        
        # B_T = self._get_magnetic_field(R_current_m);
        E_Vm = self._get_electric_field(R_current_m)
        dv_dt_0_m0s2 = get_electric_acceleration(self.q_C,self.gamma0,self.m_kg,E_Vm,self.v0_m0s,B_T=None)
        
        v_minus_half_m0s = self.v0_m0s - 1/2 * dv_dt_0_m0s2 * dt_s
        
        R_current_m = self.R0_m#np.copy(self.R0_m)
        gamma_current = self.gamma0#np.copy(self.gamma0)
        v_current_m0s = v_minus_half_m0s
        
        self.v_minus_half_vec_m0s = list([]); self.gamma_vec = list([]); self.R_vec_m = list([]); 
        self.v_minus_half_vec_m0s.append(v_current_m0s); self.gamma_vec.append(gamma_current); self.R_vec_m.append(R_current_m);

        count = 0
        R_prev_m = np.array([-np.inf,-np.inf,-np.inf])
        while self._is_in_spectrometer(R_current_m,R_prev_m):#, self.h_m, self.d_m, self.w_m, self.shield_mm, self.pinhole_dia_mm):
            count +=1
            R_star_m = R_current_m + 0.5*dt_s*v_current_m0s
            v_prev = v_current_m0s
            B_T_star = self._get_magnetic_field(R_star_m);

            E_Vm_i_star = self._get_electric_field(R_star_m)# change it
            
            v_minus_half_m0s = v_current_m0s
            gamma_i = gamma_current
            R_m_i = R_current_m 
            
            dv1_m0s= get_electric_acceleration(self.q_C,gamma_i,self.m_kg,E_Vm_i_star,v_minus_half_m0s,B_T=None) * dt_s/2
            v1_m0s = v_minus_half_m0s + dv1_m0s 
            gamma_m1 =  vel2gamma(v1_m0s)
            # R_m_m1 = R_m_i + dv1_m0s*dt_s/2

            
            
           
            
            v2_m0s  = get_magnetic_rotation(self.q_C,gamma_m1,self.m_kg,E_Vm_i_star,v1_m0s,B_T_star,dt_s)
            gamma_m2 =  vel2gamma(v2_m0s)
            
            
            v_current_m0s = v2_m0s + get_electric_acceleration(self.q_C,gamma_m2,self.m_kg,E_Vm_i_star,v2_m0s,B_T=None)*dt_s/2
            R_current_m = (R_m_i + 0.5*(v_current_m0s + v_prev)*dt_s)
            gamma_current = vel2gamma(v_current_m0s)
            
    
            self.v_minus_half_vec_m0s.append(v_current_m0s)
            self.R_vec_m.append(R_current_m)
            self.gamma_vec.append(gamma_current)
            
            R_prev_m = self.R_vec_m[-2]
            
            
        self.R_vec_m[-1], self.v_minus_half_vec_m0s[-1] ,self.gamma_vec[-1]  = self._exect_exist(self.R_vec_m, self.v_minus_half_vec_m0s)
        # self.R_vec_m[-1], self.v_minus_half_vec_m0s[-1], frac =  self._hermite_exit_crossing(self.R_vec_m[-2], self.R_vec_m[-1], self.v_minus_half_vec_m0s[-2], self.v_minus_half_vec_m0s[-1], self.dt_s, self.h_m)
        

        self.R_vec_mm = [R_m*1e3 for R_m in self.R_vec_m]
        return self.R_vec_mm,self.v_minus_half_vec_m0s,self.gamma_vec



