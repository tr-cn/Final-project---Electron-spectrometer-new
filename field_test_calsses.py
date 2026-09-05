import matplotlib.pyplot as plt
import numpy as np
import time

class Potential:
    def  __init__ (self, Ny = 2**5+1, Nx = 2**5+1, Ly_mm = 80, Lx_mm = 26, P0 = 1, start_bun = 0, end_bun = 2, gs_sweeps = 5, ie_sweeps = 4, fe_sweeps = 10, threshold = 1e-6):
        self.Ny = Ny
        self.Nx = Nx
        self.Ly_m = Ly_mm * 1e-3
        self.Lx_m = Lx_mm * 1e-3        
        
        self.dy = self.Ly_m / (Ny - 1)
        self.dx = self.Lx_m / (Nx - 1)
        self.y = np.linspace(-self.Ly_m / 2, self.Ly_m / 2, Ny)
        self.x = np.linspace(-self.Lx_m / 2, self.Lx_m / 2, Nx)

        
        self.P0 = P0
        self.start_bun = start_bun
        self.end_bun = end_bun
        self.gs_sweeps = gs_sweeps
        self.threshold = threshold
        self.ie_sweeps = ie_sweeps
        self.fe_sweeps = fe_sweeps
   
        
    def apply_boundaries (self,potential):
        # This function is overide when MagneticPotential and ElectricPotential are set
        start = self.start_bun; end = self.end_bun
        P0 = self.P0
        potential[0,:] = 0;  potential[-1,:] = 0;
        potential[0, start:end] = P0/2; potential[-1,start:end] = -P0/2;
        potential[1:,0] = 0; potential[1:,-1] = 0;
        
        return  potential
    
    
    def gauss_seidel (self, potential): 
        Ny = self.Ny; Nx = self.Nx
        
        for i in range(1, Ny - 1):
            for j in range(1, Nx - 1):
                potential[i, j] = (
                    (potential[i+1, j] + potential[i-1, j]) / self.dy**2
                    + (potential[i, j+1] + potential[i, j-1]) / self.dx**2
                ) / (2/self.dy**2 + 2/self.dx**2)
        
        return potential
    
    
    def potential_init(self):
        potential0 = np.random.rand(self.Ny, self.Nx) # [raws,colmns]
        potential0 = self.apply_boundaries (potential0)
        
        return potential0
    
    
    def potential_smoother(self, potential):
        
        for i in range(self.gs_sweeps):
            potential = self.gauss_seidel(potential)
        potential = self.apply_boundaries(potential)
        
        return potential
    
    
    def laplace_residuals (self, potential,dy,dx):
        
        r = np.zeros_like(potential)
        r[1:-1, 1:-1] =  ((potential[2:, 1:-1] - 2*potential[1:-1, 1:-1] + potential[:-2, 1:-1]) / dy**2
                         +(potential[1:-1, 2:] - 2*potential[1:-1, 1:-1] + potential[1:-1, :-2]) / dx**2)

        
        return r
    
    
    def restrict_residual (self, r):
        # r represents the residuals of finer grid
        Ny_h, Nx_h = (np.array(r.shape)+ 1)  / 2
        Ny_h = int(Ny_h); Nx_h = int(Nx_h);
        R = np.zeros([Ny_h,Nx_h])
        for I in range (1,Ny_h-1):
            for J in range (1,Nx_h-1):
                i = 2*I
                j = 2*J
                R[I,J] = 1/4 * (r [i,j] + r [i,j-1] + r [i-1,j] + r [i-1,j-1])
        r[0,:] = 0; r[-1,:] = 0; r[:,0] = 0; r[:,-1] = 0;
        
        return R
    
    
    def gauss_seidel_residuals(self, R, sweeps,dy,dx):
        # r : residuals
        # N : number of iterations
        
        Ny_h,Nx_h = R.shape
        e = np.zeros_like(R)
        for n in range (sweeps):
            for i in range (1,Ny_h-1):
                for j in range (1,Nx_h-1):
                    e[i, j] = ( (e[i + 1, j] + e[i - 1, j]) / dy**2
                            +   (e[i, j + 1] + e[i, j - 1]) / dx**2
                            -    R[i, j]) / (2 / dy**2 + 2 / dx**2)
                    
            
        e[0, :] = 0.0; e[-1, :] = 0.0; e[:, 0] = 0.0; e[:, -1] = 0.0;
        
        return e
    
    
    def restrict_potential(self, potential):
        Ny = self.Ny; Nx = self.Nx
        
        ie_sweeps = self.ie_sweeps
        fe_sweeps = self.fe_sweeps
        r = []; e = []
        dy_fine =  self.dy; dx_fine = self.dx;
        
        r_coarser_new = self.laplace_residuals(potential,dy_fine,dx_fine)
        while Ny>3 or Nx>3:
            r_finer = r_coarser_new
            r_coarser = self.restrict_residual(r_finer) 
            dy_coars = 2*dy_fine; dx_coars = 2*dx_fine;
            e_coarser = self.gauss_seidel_residuals(r_coarser, ie_sweeps,dy_coars,dx_coars)
            r_coarser_new = r_coarser - self.laplace_residuals(e_coarser,dy_coars,dx_coars)
            dy_fine = dy_coars
            dx_fine = dx_coars
            
            r.append(r_coarser_new)
            e.append(e_coarser)
            Ny, Nx = np.array(r_coarser_new.shape)
            
        
        # here the grid size is 3
        e_coarsest = self.gauss_seidel_residuals(r[-1],fe_sweeps,dy_coars,dx_coars);
        r_coarsest = r[-1] - self.laplace_residuals(e_coarsest,dy_coars,dx_coars)
        
        r[-1]=(r_coarsest)
        e[-1]=(e_coarsest)
        
        return r, e
    
    
    def prolongate_correction(self,e_coarse):
        Ny_coarse,Nx_coarse = e_coarse.shape
        Ny_fine = Ny_coarse * 2 -1
        Nx_fine = Nx_coarse * 2 -1
        e_fine = np.zeros([Ny_fine, Nx_fine])
        

        for I in range(Ny_coarse):
            for J in range(Nx_coarse):

                i = 2 * I
                j = 2 * J

                # 1. Coincident coarse/fine point
                e_fine[i, j] = e_coarse[I, J]

                # 2. Interpolation in y direction
                # f(z) = (f(2) -f(1)) / (z2-z1) * (z-z1) + f(1)  
                # The distance between z2 and z1 is 2, and between z to z1 is 1
                # f(z) = f(2) - f(1) / 2 * (1) +f(1) = (f(2) +f(1)) / 2
                if I < Ny_coarse - 1:
                    e_fine[i + 1, j] = 0.5 * (
                        e_coarse[I, J]
                        + e_coarse[I + 1, J]
                    )

                # 3. Interpolation in y direction
                if J < Nx_coarse - 1:
                    e_fine[i, j + 1] = 0.5 * (
                        e_coarse[I, J]
                        + e_coarse[I, J + 1]
                    )

                # 4. Bilinear interpolation at cell center
                if I < Ny_coarse - 1 and J < Nx_coarse - 1:
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
        potential = potential - self.prolongate_correction(e[0])

        # Restore the physical boundary conditions
        potential = self.apply_boundaries(potential)
       
        return potential
      
    
    def  laplace_func (self, potential_0):

        eps = 1;
        threshold = self.threshold
        count = 0 
        max_round = 10000
        potential = potential_0
        
        while eps > threshold and count<=max_round:
            count+=1
            
            potential = self.potential_smoother(potential)
            r, e  = self.restrict_potential(potential)
            potential = self.prolongate_potential (e, potential)
            potential = self.potential_smoother(potential)

            if count % 5 == 0:
                eps = np.max(np.abs(self.laplace_residuals(potential,self.dy,self.dx)))
                
        print (f"residual value: {eps}")
        print (f" number of iterations: {count}")
        # print(phinew)
        
        return potential
    
    
    def solve_potential (self):
        potential_0 = self.potential_init()
        self.potential   = self.laplace_func(potential_0)
        return self.potential  
            

class MagneticPotential(Potential):
    
    def __init__(self, B0_T, pole_y_start_mm, pole_y_end_mm,**kwargs):
        super().__init__(**kwargs)

        self.B0_T = B0_T
        
        self.gap_m = self.Lx_m

        # Extent of the pole faces along y
        self.pole_y_start_m = pole_y_start_mm * 1e-3
        self.pole_y_end_m = pole_y_end_mm * 1e-3

        self.mu0 =1 # 4.0 * np.pi * 1e-7
        self.delta_psi = self.B0_T * self.gap_m / self.mu0


    def apply_boundaries(self, potential):
        # The second index corresponds to y, so use self.x here.
        start = self.pole_y_start_m; end = self.pole_y_end_m;
        pole_mask = ((self.y >= start) & (self.y <= end))

        # Set all external boundaries to zero first.
        potential[0, :] = 0.0; potential[-1, :] = 0.0;
        potential[:, 0] = 0.0; potential[:, -1] = 0.0

        # Pole faces at y = ymin and y = ymax.
        potential[pole_mask,0] = +0.5 * self.delta_psi
        potential[pole_mask, -1] = -0.5 * self.delta_psi
        potential[0,1:-1] = 0; potential[-1,1:-1] = 0;
        return potential

    
class ElectricPotential(Potential):
    def __init__(self,E0_Vm, electrode_y_start_mm, electrode_y_end_mm,**kwargs):
        super().__init__(**kwargs)

        self.E0_Vm = E0_Vm
        self.gap_m = self.Lx_m
        self.electrode_y_start_m = electrode_y_start_mm * 1e-3
        self.electrode_y_end_m = electrode_y_end_mm * 1e-3
        self.voltage = E0_Vm * self.gap_m

    def apply_boundaries (self,potential):
         
         start = self.electrode_y_start_m; end = self.electrode_y_end_m
         V =  self.voltage
         electrode_mask = ((self.y >= start) & (self.y <= end))
         
         potential[0, :] = 0.0; potential[-1, :] = 0.0;
         potential[:, 0] = 0.0; potential[:, -1] = 0.0

         potential[electrode_mask,0] = +0.5 * V
         potential[electrode_mask,-1] = -0.5 * V
         
         
         potential[0,1:-1] = 0; potential[-1,1:-1] = 0;
     
         return  potential   

            
class MagneticField:
    def __init__(self, potential_object):
        self.potential_object = potential_object
        self.y = potential_object.y
        self.x = potential_object.x
        self.dy = potential_object.dy
        self.dx = potential_object.dx
        self.mu0 = 1 #4.0 * np.pi * 1e-7
        self.By = None
        self.Bx = None


    def solve_field(self):
        psi = self.potential_object.potential
        self.By = np.zeros_like(psi, dtype=float)
        self.Bx = np.zeros_like(psi, dtype=float)   
        # By = -mu0 * d(psi)/dy
        self.By[1:-1, :] = -self.mu0 * (psi[2:, :] - psi[:-2, :]) / (2.0 * self.dy)
        # Bx = -mu0 * d(psi)/dx
        self.Bx[:, 1:-1] = -self.mu0 * (psi[:, 2:] - psi[:, :-2]) / (2.0 * self.dx)

        # One-sided derivatives at the outer boundary
        self.By[0, :] = -self.mu0 * (psi[1, :] - psi[0, :]) / self.dy
        
        self.By[-1, :] = -self.mu0 * (psi[-1, :] - psi[-2, :]) / self.dy
        
        self.Bx[:, 0] = -self.mu0 * (psi[:, 1] - psi[:, 0]) / self.dx
        
        self.Bx[:, -1] = -self.mu0 * (psi[:, -1] - psi[:, -2]) / self.dx

        return self.By, self.Bx
    
    
    def show_field(self):
        
        By =  self.By; Bx = self.Bx;
        psi = self.potential_object.potential
        
        # # Coordinates of every grid point
        # Y,X = np.meshgrid(self.y,  self.x, indexing="ij")
        
        # Field magnitude
        B_magnitude = np.sqrt(By**2 + Bx**2)
        
        # Show fewer arrows, otherwise the plot is cluttered
        skip = 2
        
        plt.figure(figsize=(10, 5))
        
        plt.pcolormesh(
            self.y * 1e3,
            self.x * 1e3,
            B_magnitude.T,
            shading="auto",
            cmap="viridis",
        )
        
        plt.colorbar(label=r"$|\mathbf{B}|$ [T]")
        
        plt.streamplot(
            self.y * 1e3,
            self.x * 1e3,
            By.T,
            Bx.T,
            density=1.5,
            color="white",
            linewidth=0.8,
            arrowsize=1.2,
        )
        
        plt.xlabel("y [mm]")
        plt.ylabel("x [mm]")
        plt.title("Magnetic-field lines")
        plt.gca().set_aspect("equal")
        plt.show()
    
class ElectricField: 
    def __init__(self, potential_object):
        self.potential_object = potential_object
        self.y = potential_object.y
        self.x = potential_object.x
        self.dy = potential_object.dy
        self.dx = potential_object.dx

        self.Ey = None
        self.Ex = None

    def solve_field(self):
        phi = self.potential_object.potential

        self.Ey = np.zeros_like(phi, dtype=float)
        self.Ex = np.zeros_like(phi, dtype=float)

        # Ey = -d(phi)/dy: centered differences in the interior
        self.Ey[1:-1, :] = -(phi[2:, :] - phi[:-2, :]) / (2.0 * self.dy)

        # Ex = -d(phi)/dx: centered differences in the interior
        self.Ex[:, 1:-1] = -(phi[:, 2:] - phi[:, :-2]) / (2.0 * self.dx)

        # One-sided derivatives at the outer boundaries
        self.Ey[0, :] = -(phi[1, :] - phi[0, :]) / self.dy

        self.Ey[-1, :] = -(phi[-1, :] - phi[-2, :]) / self.dy

        self.Ex[:, 0] = -(phi[:, 1] - phi[:, 0]) / self.dx

        self.Ex[:, -1] = -(phi[:, -1] - phi[:, -2]) / self.dx

        return self.Ey, self.Ex


    def show_field(self):
        
        Ey =  self.Ey; Ex = self.Ex;
        phi = self.potential_object.potential
        E_magnitude = np.sqrt(Ey**2 + Ex**2)
        plt.figure(figsize=(10, 5))
        Y, X = np.meshgrid(self.y,  self.x, indexing="ij")
        skip = 2
        # Background: field magnitude
        plt.pcolormesh(self.y * 1e3, self.x * 1e3, E_magnitude.T, 
                       shading="auto", cmap="viridis",)
        
        plt.colorbar(label=r"$|\mathbf{E}|$ [V/m]")
        
        # Arrows: By is horizontal, Bx is vertical
        plt.streamplot( self.y * 1e3, self.x * 1e3, Ey.T, Ex.T, density=1.5,
            color="white", linewidth=0.8, arrowsize=1.2)
        
        plt.xlabel("y [mm]")
        plt.ylabel("x [mm]")
        plt.title("Electric field magnitude and direction")
        plt.gca().set_aspect("equal")
        plt.show()
        

if __name__ == "__main__":
    plt.close('all')
    Ny = 2**5+1; Nx = 2**5+1;
    V0 = 10
    B0 = 0.5
    
    
    phi = ElectricPotential(Ny=2**5+1, Nx=2**5+1, Ly_mm=12.5, Lx_mm=12.5, E0_Vm=10,
                            electrode_y_start_mm=-10, electrode_y_end_mm= 10)
    
    phi.solve_potential()
    Ele = ElectricField(phi)
    Ele.solve_field()
    Ele.show_field()
    
    
    if True:
        psi = MagneticPotential(Ny=2**5+1, Nx=2**6+1, Ly_mm=80, Lx_mm=12.5, B0_T=10,
                                pole_y_start_mm=-40.0, pole_y_end_mm=20.0)
        
        psi.solve_potential()
        mag = MagneticField(psi)
        mag.solve_field()
        mag.show_field()
        


