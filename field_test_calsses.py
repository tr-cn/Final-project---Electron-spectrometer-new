import matplotlib.pyplot as plt
import numpy as np
import time

class Potential:
    def  __init__ (self, Nx = 2**5+1, Ny = 2**5+1, P0 = 1, start_bun = 0, end_bun = 2, gs_swips = 5, ie_swips = 4, fe_swips = 10, treshhold = 1e-6):
        self.Nx = Nx
        self.Ny = Ny
        self.P0 = P0
        self.start_bun = start_bun
        self.end_bun = end_bun
        self.gs_swips = gs_swips
        self.treshhold = treshhold
        self.ie_swips = ie_swips
        self.fe_swips = fe_swips
        
    def apply_boundries (self,potential):
        start = self.start_bun; end = self.end_bun
        P0 = self.P0
        potential[0,:] = 0;  potential[-1,:] = 0;
        potential[0, start:end] = P0/2; potential[-1,start:end] = -P0/2;
        potential[1:,0] = 0; potential[1:,-1] = 0;
        return  potential
    
    def gauss_seidel (self, potential): 
        Nx = self.Nx; Ny = self.Ny
        
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                potential[i,j] = 1/4 * (potential[i+1,j] + potential[i-1,j] + potential[i,j+1] + potential[i,j-1])
        potential = self.apply_boundries (potential)
        
        return potential
    
    def potential_init(self):
        potential0 = np.random.rand(self.Nx, self.Ny) # [raws,colmns]
        potential0 = self.apply_boundries (potential0)
        return potential0
    
    def potential_smother(self, potential):
        
        for i in range(self.gs_swips):
            potential = self.gauss_seidel(potential)
        potential = self.apply_boundries(potential)
        return potential
    
    def laplace_residuals (self, potential):
        r = np.zeros_like(potential)
        r[1:-1, 1:-1] = (potential[2:, 1:-1] + potential[:-2, 1:-1] + potential[1:-1, 2:] + potential[1:-1, :-2] - 4.0 * potential[1:-1, 1:-1])
        
        return r
    
    
    def restrict_residual (self, r):
        # r represents the residuals of finer grid
        Nx_h, Ny_h = (np.array(r.shape)+ 1)  / 2
        Nx_h = int(Nx_h); Ny_h = int(Ny_h);
        R = np.zeros([Nx_h,Ny_h])
        for I in range (1,Nx_h-1):
            for J in range (1,Ny_h-1):
                i = 2*I
                j = 2*J
                R[I,J] = 1/4 * (r [i,j] + r [i,j-1] + r [i-1,j] + r [i-1,j-1])
        r[0,:] = 0; r[-1,:] = 0; r[:,0] = 0; r[:,-1] = 0;
        
        return R
    
    def gauss_seidel_residuals(self, R, swips):
        # r : residuals
        # N : number of iterations
        Nx_h,Ny_h = R.shape
        
        e = np.zeros_like(R)
        for n in range (swips):
            for i in range (1,Nx_h-1):
                for j in range (1,Ny_h-1):
                    e[i,j] = 1/4 * (e[i+1,j] + e[i-1,j] + e[i,j+1] + e[i,j-1] + R[i,j])
            
        e[0, :] = 0.0; e[-1, :] = 0.0; e[:, 0] = 0.0; e[:, -1] = 0.0;
        
        return e
    
    def restrict_potential(self, potential):
        Nx = self.Nx; Ny = self.Ny
        
        ie_swips = self.ie_swips
        fe_swips = self.fe_swips
        r = []; e = []
            
        r.append(self.laplace_residuals(potential))
        while Nx>3 or Ny>3:
            r_finer = r[-1]
            r_coarser = self.restrict_residual(r_finer)
            e_coarser = self.gauss_seidel_residuals(r_coarser, ie_swips)
            r_coarser_new = r_coarser - self.laplace_residuals(e_coarser)
           
            r.append(r_coarser_new)
            e.append(e_coarser)
            Nx, Ny = np.array(r_coarser_new.shape)
        
        # here the grid size is 3
        e_coarsest = self.gauss_seidel_residuals(r[-1],fe_swips);
        r_coarsest = r_coarser - self.laplace_residuals(e_coarsest)
        
        r[-1]=(r_coarsest)
        e[-1]=(e_coarsest)
        return r, e
    
    def prolongate_correction(self,e_coarse):
        Nx_coarse,Ny_coarse = e_coarse.shape
        Nx_fine = Nx_coarse * 2 -1
        Ny_fine = Ny_coarse * 2 -1
        e_fine = np.zeros([Nx_fine, Ny_fine])
        

        for I in range(Nx_coarse):
            for J in range(Ny_coarse):

                i = 2 * I
                j = 2 * J

                # 1. Coincident coarse/fine point
                e_fine[i, j] = e_coarse[I, J]

                # 2. Interpolation in x direction
                # f(x) = (f(2) -f(1)) / (x2-x1) * (x-x1) + f(1)  
                # The distance between x2 and x1 is 2, and between x to x1 is 1
                # f(x) = f(2) - f(1) / 2 * (1) +f(1) = (f(2) +f(1)) / 2
                if I < Nx_coarse - 1:
                    e_fine[i + 1, j] = 0.5 * (
                        e_coarse[I, J]
                        + e_coarse[I + 1, J]
                    )

                # 3. Interpolation in y direction
                if J < Ny_coarse - 1:
                    e_fine[i, j + 1] = 0.5 * (
                        e_coarse[I, J]
                        + e_coarse[I, J + 1]
                    )

                # 4. Bilinear interpolation at cell center
                if I < Nx_coarse - 1 and J < Ny_coarse - 1:
                    e_fine[i + 1, j + 1] = 0.25 * (
                        e_coarse[I, J]
                        + e_coarse[I + 1, J]
                        + e_coarse[I, J + 1]
                        + e_coarse[I + 1, J + 1]
                    )

        # Correction has zero Dirichlet boundary values
        e_fine[0, :] = 0.0; e_fine[-1, :] = 0.0; e_fine[:, 0] = 0.0; e_fine[:, -1] = 0.0

        return e_fine

    def prolongate_potential(self,e, potential):
        # Go from the coarsest correction upward
        for level in range(len(e) - 2, -1, -1):
            finer_correction = self.prolongate_correction(e[level + 1])

            # Add the correction to the current level
            e[level] = e[level] + finer_correction

        # Transfer the finest correction to potential
        potential = potential + self.prolongate_correction(e[0])

        # Restore the physical boundary conditions
        potential = self.apply_boundries(potential)
       

        return potential
        
    def  laplace_func (self, potential_0):
        Nx = self.Nx; Ny = self.Ny;
        gs_swips = self.gs_swips # gauss-seidle swips befor start mutigrid
        eps = 1;
        treshhold = self.treshhold
        count = 0 
        max_round = 10000
        potential = potential_0
        
        while eps > treshhold and count<=max_round:
            potential_prev = np.copy(potential)
            count+=1
            
            potential = self.potential_smother(potential)
            r, e  = self.restrict_potential(potential)
            potential = self.prolongate_potential (e, potential)
            potential = self.potential_smother(potential)
          
            if count % 5 == 0:
                eps = np.max(np.abs(self.laplace_residuals(potential)))
                
        print (f"residual value: {eps}")
        print (f" number of iterations: {count}")
        # print(phinew)
        
        return potential
    
    def get_potential (self):
        potential_0 = self.potential_init()
        potential   = self.laplace_func(potential_0)
        
        return potential

# class field:
    
    


if __name__ == "__main__":
    plt.close('all')
    Nx = 2**5+1; Ny = 2**5+1;
    V0 = 10
    B0 = 0.5
    
    
    phi = Potential (Nx = 2**5+1, Ny = 2**5+1, P0 = 1, start_bun = 2**3, end_bun = 2**4)
    psi = Potential (Nx = 2**5+1, Ny = 2**5+1, P0 = 1, start_bun = 0, end_bun = 2**4)
    plt.figure()
    phi_p = phi.get_potential()
    plt.imshow(phi_p,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    plt.show()
    
    plt.figure()
    psi_p = psi.get_potential()
    plt.imshow(psi_p,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    plt.show()
