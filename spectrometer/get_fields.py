import numpy as np


class Magenetic_field_Analitic:
    def __init__(self, B0_T, width_mm=12.5, depth_mm=50.8, k=0, fringe = 1, sharp_edge=1, yoke=1):
        self.B0_T = B0_T
        self.k = k # Gradient along Z direction
        self.fringe = fringe
        self.yoke = yoke
        
        self.width_m = width_mm*1e-3
        self.depth_m = depth_mm*1e-3
        if sharp_edge:
            self.a0 = 0.3835
            self.a1 = 2.388
            self.a2 = -0.8171
            self.a3 = 0.2 
        else:
            self.a0 = 0.531
            self.a1 = 2.341
            self.a2 = 0.7799
            self.a3 = 0.110
            
        
    def _enge_func(self, y_complex):
        S = self.a0 + self.a1 * y_complex + self.a2 * (y_complex**2) + self.a3 * (y_complex**3)
        return 1 / (1 + np.exp(S))
    
    def _get_magnetic_field(self,R):
        if self.fringe == 0:
            return np.array([self.B0_T,0,0])
        
        
        x = R[0]; y = R[1]; z = R[2];
        
        y_complex = y + 1j*x
        y_comlex_norm = -y_complex/self.width_m
        F_in = self._enge_func(y_comlex_norm)
        
        if self.yoke:
            F = F_in
        else:
            y_complex_norm_out = (y_complex - self.depth_m)/self.width_m
            F_out = self._enge_func(y_complex_norm_out)
            F = F_in*F_out
        
        Bx = self.B0_T * np.real(F) * (1 - self.k * z)
        By = self.B0_T * np.imag(F)
        Bz = self.B0_T * self.k * x
        
        return np.array([Bx,By,Bz])