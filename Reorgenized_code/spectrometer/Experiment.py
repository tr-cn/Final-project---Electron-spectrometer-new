import numpy as np

class Experiment():
    def __init__(self, Spec_class, q_eng_MeV, m_kg, q_C,Bx_T, Ey_Vm):
        self.width_mm = Spec_class.width_mm
        self.height_mm = Spec_class.height_mm
        self.depth_mm = Spec_class.depth_mm
        self.shield_mm = Spec_class.shield_mm
        self.m_kg = m_kg
        self.q_C = q_C
        self.Bx_T = Bx_T
        self.Ey_Vm = Ey_Vm
        
        
    
     def _time_grid():
         
         
        
        
    def _experiment_bonds(self,R0_mm, grid, fringe = 1):
        Nx = grid[0]; Ny = grid[1]; Nz = grid[2];
        self.Lx_mm = np.linspace(-self.width_mm/2,self.width_mm/2,Nx)
        self.Lz_mm = np.linspace(-self.height_mm/2,self.height_mm/2,Nz)
        if self.yoke == 1:
            self.Ly_mm = np.linsapce(-abs(self.shield_mm),self.depth_mm,Ny)
        elif self.yoke == 0:
            self.Ly_mm = np.linsapce(-abs(self.shield_mm),self.depth_mm + abs(self.shield_mm),Ny)
            
        return self.Lx_mm, self.Ly_mm, self.Lz_mm
            