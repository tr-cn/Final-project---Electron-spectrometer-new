import numpy as np


def MeV2m0s(q_eng_MeV, m_kg):
    q_eng_J = abs(q_eng_MeV)*1e6 * 1.602*1e-19
    c_m0s = 299792458
    v0_m0s = c_m0s*np.sqrt ( 1 - ( m_kg*c_m0s**2 / (q_eng_J + m_kg*c_m0s**2) )**2 )# [m/s]
    return v0_m0s

def vel2gamma (v_m0s):
    c_m0s = 299792458
    v0_mag_m0s = np.linalg.norm(v_m0s)
    return 1 / np.sqrt ( 1 - (v0_mag_m0s/c_m0s)**2 )




class Experiment():
    def __init__(self, height_mm, width_mm, depth_mm, shield_mm, yoke,
                 q_eng_MeV, m_kg, q_C,
                 Bx0_T, Ey0_Vm,
                 N_steps = 150, CFL=0.1):
        
        self.width_mm = width_mm
        self.height_mm = height_mm
        self.depth_mm = depth_mm
        self.shield_mm = shield_mm
        self.yoke = yoke
        
        self.q_eng_MeV = q_eng_MeV
        self.m_kg = m_kg
        self.q_C = q_C
        self.Bx0_T = Bx0_T
        self.Ey0_Vm = Ey0_Vm
        
        self.CFL = CFL
        self.N_steps = N_steps
        
        
        
    def _get_T_cyclotorn(self):
        self.v_m0s = MeV2m0s (self.q_eng_MeV,self.m_kg)
        self.gamma = vel2gamma (self.v_m0s)
        self.T_cyclotron =  ( abs(self.q_C)*np.linalg.norm(self.Bx0_T)/(np.pi*self.gamma*self.m_kg) )**-1   
        return  self.T_cyclotron

    def _get_dt(self):
        
        v0_m0s = MeV2m0s(self.q_eng_MeV, self.m_kg)
        T_cyclotron = self._get_T_cyclotorn()
        self.dt =T_cyclotron/self.N_steps

        return self.dt
         
    
    def _get_dr(self):
        #  I call it that way, althogh it is not propagation of a wave
        self.CFL = CFL
        self.dr = self.dt*self.CFL;
        return self.dr 
    
    def _get_experiment_range(self):
        self.x_start_mm = -self.width_mm/2
        self.x_end_mm   =  self.width_mm/2
        
        self.z_start_mm = -self.height_mm/2
        self.z_end_mm   =  self.height_mm/2
        
        self.y_start_mm = -abs(self.shield_mm)
        if self.yoke == 1:
            self.y_end_mm = self.depth_mm
                 
        elif self.yoke == 0:
            self.y_end_mm + abs(self.shield_mm)
    
        self.exp_range_X_mm = np.array([self.x_start_mm, self.x_end_mm])
        self.exp_range_Y_mm = np.array([self.y_start_mm, self.y_end_mm])
        self.exp_range_Z_mm = np.array([self.z_start_mm, self.z_end_mm])
        return self.exp_range_X_mm, self.exp_range_Y_mm, self.exp_range_Z_mm
    
    
    def _get_diffs(self):
        exp_range_X_mm, exp_range_Y_mm, exp_range_Z_mm = self._get_experiment_range()
        
        Dx_mm = np.diff(exp_range_X_mm)
        Dy_mm = np.diff(exp_range_Y_mm)
        Dz_mm = np.diff(exp_range_Z_mm)
        
        return Dx_mm, Dy_mm, Dz_mm
        
    def _get_dxdydz(self):
        Dx_mm, Dy_mm, Dz_mm = self._get_diffs()
        
        norm = np.sqrt(Dx_mm**2 + Dy_mm**2 + Dz_mm**2)
        self.dx_mm  = self.dr * Dx_mm/norm
        self.dy_mm  = self.dr * Dy_mm/norm
        self.dz_mm  = self.dr * Dz_mm/norm
        
        
        
        return self.dx_mm, self.dy_mm, self.dz_mm
        
    
    def _get_needed_grid_size(self):
        Dx_mm, Dy_mm, Dz_mm = self._get_diffs()
        norm = np.sqrt(Dx_mm**2 + Dy_mm**2 + Dz_mm**2)
        N = self.N_steps/self.CFL
        self.Nx = int(np.round(N*Dx_mm/norm)); 
        self.Ny = int(np.round(N*Dy_mm/norm));
        self.Nz = int(np.round(N*Dz_mm/norm));
        return self.Nx, self.Ny, self.Nz
                          
                           
    def _get_range_vecs(self):
        self.Lx_mm = np.linspace(self.x_start_mm, self.x_end_mm,self.Nx)
        self.Ly_mm = np.linspace(self.y_start_mm, self.y_end_mm,self.Ny)
        self.Lz_mm = np.linspace(self.z_start_mm, self.z_end_mm,self.Nz)
        return  self.Lx_mm,  self.Ly_mm,  self.Lz_mm
    
    
    
    def _evaluate_exp_paramas(self):
        self._get_T_cyclotorn()
        self._get_dt()
        self._get_dr()
        self._get_experiment_range()
        self._get_diffs()
        self._get_dxdydz()
        self._get_needed_grid_size()
        self._get_range_vecs()
        
        
        
        
    # def _experiment_bondriess(self,R0_mm, grid, fringe = 1):
    #     Nx = grid[0]; Ny = grid[1]; Nz = grid[2];
    #     self.Lx_mm = np.linspace(-self.width_mm/2,self.width_mm/2,Nx)
    #     self.Lz_mm = np.linspace(-self.height_mm/2,self.height_mm/2,Nz)
    #     if self.yoke == 1:
    #         self.Ly_mm = np.linsapce(-abs(self.shield_mm),self.depth_mm,Ny)
    #     elif self.yoke == 0:
    #         self.Ly_mm = np.linsapce(-abs(self.shield_mm),self.depth_mm + abs(self.shield_mm),Ny)
            
    #     return self.Lx_mm, self.Ly_mm, self.Lz_mm
if __name__ == "__main__":
    me_kg = 9.109*1e-31
    Bx0_T = 0.5
    Ey0_Vm = 0
    E_V0m = np.array([0,0,0])
    R0_mm  = np.array([0,-12.5,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    steps = 150
    yoke = 1
    shield_mm = 12.5
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.1
    N_steps = 150
    experiment = Experiment(height_mm = height_mm, width_mm = width_mm, depth_mm = depth_mm,  q_eng_MeV = e_eng_MeV, m_kg = me_kg, 
                            q_C = e_C ,Bx0_T = Bx0_T, Ey0_Vm = Ey0_Vm, yoke = yoke, shield_mm = shield_mm ,CFL =CFL, N_steps = N_steps)
       
    experiment._evaluate_exp_paramas()