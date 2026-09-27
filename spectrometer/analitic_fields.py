import numpy as np


# def 


class Magenetic_field_Analitic:
    def __init__(self, Bx0_T, width_mm=12.5, depth_mm=50.8, k=0, fringe = 1, sharp_edge=1, yoke=1):
        self.Bx0_T = Bx0_T
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
            return np.array([self.Bx0_T,0,0])
        
        
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
        
        Bx = self.Bx0_T * np.real(F) * (1 - self.k * z)
        By = self.Bx0_T * np.imag(F)
        Bz = self.Bx0_T * self.k * x
        
        return np.array([Bx,By,Bz])
        
    
    
class Electric_field_Analitic:
    def __init__(self, Ex0_Vm, width_mm=12.5, depth_mm=50.8, k=0, fringe=1, sharp_edge=1, yoke=1):
        self.Ex0_Vm = Ex0_Vm
        self.k = k
        self.fringe = fringe
        self.yoke = yoke
        self.width_m = width_mm*1e-3
        self.depth_m = depth_mm*1e-3
        if sharp_edge:
            self.a0, self.a1, self.a2, self.a3 = 0.3835, 2.388, -0.8171, 0.2
        else:
            self.a0, self.a1, self.a2, self.a3 = 0.531, 2.341, 0.7799, 0.110

    def _enge_func(self, y_complex):
        S = self.a0 + self.a1*y_complex + self.a2*(y_complex**2) + self.a3*(y_complex**3)
        return 1/(1+np.exp(S))

    def _get_electric_field(self, R):
        if self.fringe == 0:
            return np.array([0, self.Ex0_Vm, 0])

        x, y, z = R[0], R[1], R[2]
        y_complex = y + 1j*x
        y_complex_norm = -y_complex/self.width_m
        F_in = self._enge_func(y_complex_norm)

        if self.yoke:
            F = F_in
        else:
            y_complex_norm_out = (y_complex - self.depth_m)/self.width_m
            F_out = self._enge_func(y_complex_norm_out)
            F = F_in*F_out

        Ex = self.Ex0_Vm * np.real(F)
        Ey = self.Ex0_Vm * np.imag(F) * (1 - self.k*z)
        Ez = self.Ex0_Vm * self.k * x
        return np.array([Ex, Ey, Ez])
    
    

# class Electric_field:
#     def __init__(self, V, width_mm=12.5, depth_mm=50.8, mesh_yz = [128,128], k = 0, fringe = 1, yoke=1):
#         self.V = V
#         self.mesh_yz = mesh_yz;
#         self.k = k # Gradient along Z direction
#         self.fringe = fringe
#         self.yoke = yoke
#         self.width_m = width_mm*1e-3
#         self.depth_m = depth_mm*1e-3
        
#     def phi_init (self):
#         Ny = self.mesh_yz[0]        
#         Nz =self.mesh_yz[1]
#         V = self.V
#         phi = np.linspace(V, 0, Nz).reshape(-1, 1) * np.ones((1, Ny))
#         return phi
    
        
#     def ghost_cell (self):
#         # Derichle 
#         phi = self.phi
#         phi_g =  np.zeros(np.array(phi.shape)+2)
#         phi_g[1:-1,1:-1] = phi.copy()
#         self.phi_g = phi_g
#         return phi_g
        
    
#     def residual (self,phi_prev,phi_new):
#         eps_max = np.max(np.max(np.abs(phi_new - phi_prev)))
#         # print(np.max(phi_prev))
#         # print(np.min((phi_new - phi_prev)**2))
#         return eps_max
    
    
#     def laplace_func (self):
#         phi = self.phi_init();
#         phi_new = np.zeros_like(phi)
#         Ny = self.mesh_yz[0]        
#         Nz =self.mesh_yz[1]
#         eps = np.inf
#         count = 0;
#         while eps>1-2:
#             if ~np.mod(count,5):
#                 for i in range(1,Ny-1):
#                    for j in range(1,Nz-1):
#                        phi_new[i,j] = 1/4 * (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])  # No charges within the capacitor
#             # else:
#             #     while Nx or Ny == 8
#             #     phi_re
               
            
               
        
        
        
        
    
    
    
    
    
    
    
    
    # def _enge_func(self, y_complex):
    #     S = self.a0 + self.a1 * y_complex + self.a2 * (y_complex**2) + self.a3 * (y_complex**3)
    #     return 1 / (1 + np.exp(S))
    
    # def _get_magnetic_field(self,R):
    #     if self.fringe == 0:
    #         return np.array([self.Bx0_T,0,0])
        
        
    #     x = R[0]; y = R[1]; z = R[2];
        
    #     y_complex = y + 1j*x
    #     y_comlex_norm = -y_complex/self.width_m
    #     F_in = self._enge_func(y_comlex_norm)
        
    #     if self.yoke:
    #         F = F_in
    #     else:
    #         y_complex_norm_out = (y_complex - self.depth_m)/self.width_m
    #         F_out = self._enge_func(y_complex_norm_out)
    #         F = F_in*F_out
        
    #     Bx = self.Bx0_T * np.real(F) * (1 - self.k * z)
    #     By = self.Bx0_T * np.imag(F)
    #     Bz = self.Bx0_T * self.k * x
        
    #     return np.array([Bx,By,Bz])

        
        
