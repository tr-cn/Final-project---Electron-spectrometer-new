import matplotlib.pyplot as plt
import numpy as np
import time

class Potential:
    def  __init__ (self, Nx = 2**5+1, Ny = 2**5+1, Lx_mm = 80, Ly_mm = 26, P0 = 1, start_bun = 0, end_bun = 2, gs_sweeps = 5, ie_sweeps = 4, fe_sweeps = 10, threshold = 1e-6):
        self.Nx = Nx
        self.Ny = Ny
        self.Lx_m = Lx_mm * 1e-3
        self.Ly_m = Ly_mm * 1e-3        
        
        self.dx = self.Lx_m / (Nx - 1)
        self.dy = self.Ly_m / (Ny - 1)
        self.x = np.linspace(-self.Lx_m / 2, self.Lx_m / 2, Nx)
        self.y = np.linspace(-self.Ly_m / 2, self.Ly_m / 2, Ny)

        
        self.P0 = P0
        self.start_bun = start_bun
        self.end_bun = end_bun
        self.gs_sweeps = gs_sweeps
        self.threshold = threshold
        self.ie_sweeps = ie_sweeps
        self.fe_sweeps = fe_sweeps
   
        
    def apply_boundaries (self,potential):
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
        potential = self.apply_boundaries (potential)
        
        return potential
    
    
    def potential_init(self):
        potential0 = np.random.rand(self.Nx, self.Ny) # [raws,colmns]
        potential0 = self.apply_boundaries (potential0)
        
        return potential0
    
    
    def potential_smoother(self, potential):
        
        for i in range(self.gs_sweeps):
            potential = self.gauss_seidel(potential)
        potential = self.apply_boundaries(potential)
        
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
    
    
    def gauss_seidel_residuals(self, R, sweeps):
        # r : residuals
        # N : number of iterations
        Nx_h,Ny_h = R.shape
        
        e = np.zeros_like(R)
        for n in range (sweeps):
            for i in range (1,Nx_h-1):
                for j in range (1,Ny_h-1):
                    e[i,j] = 1/4 * (e[i+1,j] + e[i-1,j] + e[i,j+1] + e[i,j-1] + R[i,j])
            
        e[0, :] = 0.0; e[-1, :] = 0.0; e[:, 0] = 0.0; e[:, -1] = 0.0;
        
        return e
    
    
    def restrict_potential(self, potential):
        Nx = self.Nx; Ny = self.Ny
        
        ie_sweeps = self.ie_sweeps
        fe_sweeps = self.fe_sweeps
        r = []; e = []
            
        r.append(self.laplace_residuals(potential))
        while Nx>3 or Ny>3:
            r_finer = r[-1]
            r_coarser = self.restrict_residual(r_finer)
            e_coarser = self.gauss_seidel_residuals(r_coarser, ie_sweeps)
            r_coarser_new = r_coarser - self.laplace_residuals(e_coarser)
           
            r.append(r_coarser_new)
            e.append(e_coarser)
            Nx, Ny = np.array(r_coarser_new.shape)
        
        # here the grid size is 3
        e_coarsest = self.gauss_seidel_residuals(r[-1],fe_sweeps);
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
        potential = self.apply_boundaries(potential)
       
        return potential
      
    
    def  laplace_func (self, potential_0):
        Nx = self.Nx; Ny = self.Ny;
        gs_sweeps = self.gs_sweeps # gauss-seidle sweeps befor start mutigrid
        eps = 1;
        threshold = self.threshold
        count = 0 
        max_round = 10000
        potential = potential_0
        
        while eps > threshold and count<=max_round:
            potential_prev = np.copy(potential)
            count+=1
            
            potential = self.potential_smoother(potential)
            r, e  = self.restrict_potential(potential)
            potential = self.prolongate_potential (e, potential)
            potential = self.potential_smoother(potential)
          
            if count % 5 == 0:
                eps = np.max(np.abs(self.laplace_residuals(potential)))
                
        print (f"residual value: {eps}")
        print (f" number of iterations: {count}")
        # print(phinew)
        
        return potential
    
    
    def solve_potential (self):
        potential_0 = self.potential_init()
        self.potential   = self.laplace_func(potential_0)
        return self.potential


# class MagneticPotential(Potential):
#     def __init__(self, B0_T, gap_mm, pole_x_start_mm, pole_x_end_mm, **kwargs):
#         super().__init__(**kwargs)
        
#         self.B0_T = B0_T
#         self.gap_m = gap_m = gap_mm * 1e-3
#         self.pole_x_start_m = pole_x_start_mm * 1e-3
#         self.pole_x_end_m = pole_x_end_mm * 1e-3
#         mu0 = 4.0 * np.pi * 1e-7
#         self.delta_psi = B0_T * gap_m / mu0
        
#     def apply_boundaries (self,potential):
        
#         start = self.pole_x_start_m; end = self.pole_x_end_m
#         d_psi =  self.delta_psi
           
#         pole_mask = ((self.x >= start) & (self.x <= end))
#         potential[0,:] = 0;  potential[-1,:] = 0;
#         # potential[0, start:end] = V/2; potential[-1,start:end] = -V/2;
#         potential[0 , pole_mask] = +0.5 * d_psi
#         potential[-1, pole_mask] = -0.5 * d_psi
           
#         potential[1:,0] = 0; potential[1:,-1] = 0;
        
#         return  potential   
            

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
        # The second index corresponds to y, so use self.y here.
        pole_mask = ((self.y >= self.pole_y_start_m) & (self.y <= self.pole_y_end_m))

        # Set all external boundaries to zero first.
        potential[0, :] = 0.0
        potential[-1, :] = 0.0
        potential[:, 0] = 0.0
        potential[:, -1] = 0.0

        # Pole faces at x = xmin and x = xmax.
        potential[0, pole_mask] = +0.5 * self.delta_psi
        potential[-1, pole_mask] = -0.5 * self.delta_psi

        return potential


    
class ElectricPotential(Potential):
    def __init__(self,E0, gap, electrode_x_start_mm, electrode_x_end_mm,**kwargs):
        super().__init__(**kwargs)

        self.E0 = E0
        self.plate_gap = gap
        self.electrode_x_start_m = electrode_x_start_mm * 1e-3
        self.electrode_x_end_m = electrode_x_end_mm * 1e-3
        self.voltage = E0 * gap

    def apply_boundaries (self,potential):
         
         start = self.electrode_x_start_m; end = self.electrode_x_end_m
         V =  self.voltage
         electrode_mask = ((self.x >= start) & (self.x <= end))
         potential[0,:] = 0;  potential[-1,:] = 0;
         # potential[0, start:end] = V/2; potential[-1,start:end] = -V/2;
         potential[0 , electrode_mask] = +0.5 * V
         potential[-1, electrode_mask] = -0.5 * V
         
         
         potential[1:,0] = 0; potential[1:,-1] = 0;
     
         return  potential   
            


class MagneticField:
    def __init__(self, potential_object):
        self.potential_object = potential_object
        self.x = potential_object.x
        self.y = potential_object.y
        self.dx = potential_object.dx
        self.dy = potential_object.dy
        self.mu0 = 1 #4.0 * np.pi * 1e-7
        self.Bx = None
        self.By = None


    def solve_field(self):
        psi = self.potential_object.potential
        self.Bx = np.zeros_like(psi, dtype=float)
        self.By = np.zeros_like(psi, dtype=float)   
        # Bx = -mu0 * d(psi)/dx
        self.Bx[1:-1, :] = -self.mu0 * (psi[2:, :] - psi[:-2, :]) / (2.0 * self.dx)
        # By = -mu0 * d(psi)/dy
        self.By[:, 1:-1] = -self.mu0 * (psi[:, 2:] - psi[:, :-2]) / (2.0 * self.dy)

        # One-sided derivatives at the outer boundary
        self.Bx[0, :] = -self.mu0 * (psi[1, :] - psi[0, :]) / self.dx
        
        self.Bx[-1, :] = -self.mu0 * (psi[-1, :] - psi[-2, :]) / self.dx
        
        self.By[:, 0] = -self.mu0 * (psi[:, 1] - psi[:, 0]) / self.dy
        
        self.By[:, -1] = -self.mu0 * (psi[:, -1] - psi[:, -2]) / self.dy

        return self.Bx, self.By
    
    
    def show_field(self):
            
    
        Bx =  self.Bx; By = self.By;
        psi = self.potential_object.potential
        
        # Coordinates of every grid point
        X, Y = np.meshgrid(self.x,  self.y, indexing="ij")
        
        # Field magnitude
        B_magnitude = np.sqrt(Bx**2 + By**2)
        
        # Show fewer arrows, otherwise the plot is cluttered
        skip = 2
        
        plt.figure(figsize=(10, 5))
        
        # Background: field magnitude
        plt.pcolormesh(
            self.x * 1e3,
            self.y * 1e3,
            B_magnitude.T,
            shading="auto",
            cmap="viridis",
        )
        
        plt.colorbar(label=r"$|\mathbf{B}|$ [T]")
        
        # Arrows: Bx is horizontal, By is vertical
        plt.quiver(
            X[::skip, ::skip] * 1e3,
            Y[::skip, ::skip] * 1e3,
            Bx[::skip, ::skip],
            By[::skip, ::skip],
            color="white",
            angles="xy",
            scale_units="xy",
            scale=None,
            width=0.003,
        )
        
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
        plt.title("Magnetic field magnitude and direction")
        plt.gca().set_aspect("equal")
        plt.show()
        
        
        plt.figure(figsize=(10, 5))
        
        plt.pcolormesh(
            self.x * 1e3,
            self.y * 1e3,
            B_magnitude.T,
            shading="auto",
            cmap="viridis",
        )
        
        plt.colorbar(label=r"$|\mathbf{B}|$ [T]")
        
        plt.streamplot(
            self.x * 1e3,
            self.y * 1e3,
            Bx.T,
            By.T,
            density=1.5,
            color="white",
            linewidth=0.8,
            arrowsize=1.2,
        )
        
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
        plt.title("Magnetic-field lines")
        plt.gca().set_aspect("equal")
        plt.show()
                
            
        

class ElectricField: 
    def __init__(self, potential_object):
        self.potential_object = potential_object
        self.x = potential_object.x
        self.y = potential_object.y
        self.dx = potential_object.dx
        self.dy = potential_object.dy

        self.Ex = None
        self.Ey = None

    def solve_field(self):
        phi = self.potential_object.potential

        self.Ex = np.zeros_like(phi, dtype=float)
        self.Ey = np.zeros_like(phi, dtype=float)

        # Ex = -d(phi)/dx: centered differences in the interior
        self.Ex[1:-1, :] = -(phi[2:, :] - phi[:-2, :]) / (2.0 * self.dx)

        # Ey = -d(phi)/dy: centered differences in the interior
        self.Ey[:, 1:-1] = -(phi[:, 2:] - phi[:, :-2]) / (2.0 * self.dy)

        # One-sided derivatives at the outer boundaries
        self.Ex[0, :] = -(phi[1, :] - phi[0, :]) / self.dx

        self.Ex[-1, :] = -(phi[-1, :] - phi[-2, :]) / self.dx

        self.Ey[:, 0] = -(phi[:, 1] - phi[:, 0]) / self.dy

        self.Ey[:, -1] = -(phi[:, -1] - phi[:, -2]) / self.dy

        return self.Ex, self.Ey


    def show_field(self):
        
        X, Y = np.meshgrid(self.x, self.y, indexing="ij")
        
        skip = 2
        
        plt.figure(figsize=(8, 4))
        
        plt.quiver(
        X[::skip, ::skip],
        Y[::skip, ::skip],
        self.Ex[::skip, ::skip],
        self.Ey[::skip, ::skip],
        angles="xy",
        scale_units="xy",
        scale=None,
        )
        
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.title("Electric field vector")
        plt.gca().set_aspect("equal")
        plt.show()


if __name__ == "__main__":
    plt.close('all')
    Nx = 2**5+1; Ny = 2**5+1;
    V0 = 10
    B0 = 0.5
    
       
    
    # psi = MagneticPotential(
    # Nx=65,
    # Ny=65,
    # Lx=0.10,
    # Ly=0.06,
    # B0=0.5,
    # gap=0.01,
    # pole_x_start=-0.02,
    # pole_x_end=0.02,
    # )


   
#     phi = ElectricPotential (Nx = 2**5+1, Ny = 2**5+1, Lx_mm = 100, Ly_mm = 30, E0 = 1, gap = 12.5, electrode_x_start = -25, electrode_x_end = 25); phi.solve_potential()
#     plt.imshow(
#     phi.potential,
#     extent=[phi.x[0]*1e3, phi.x[-1]*1e3, phi.y[0]*1e3, phi.y[-1]*1e3],
#     origin="lower",
#     aspect="auto",
# )
# plt.colorbar()
    # plt.imshow(
    # phi.potential.T,
    # extent=[psi.x[0]*1e3, psi.x[-1]*1e3, psi.y[0]*1e3, psi.y[-1]*1e3],
    # origin="lower",
    # aspect="auto",
# )
# plt.colorbar()
    # plt.show()
   
    
    # psi = MagneticPotential(Nx = 2**5+1, Ny = 2**5+1, B0 = 1, gap = 0.05, pole_start = 2**3, pole_end = 2**4); psi.solve_potential()
    # psi.solve_potential()
    
    # magnetic_field = MagneticField(psi)
    # Bx, By = magnetic_field.solve_field()
    
    
    
    # B_magnitude = np.sqrt(Bx**2 + By**2)

    # X, Y = np.meshgrid(psi.x, psi.y, indexing="ij")
    # skip = 2
    
    # plt.figure(figsize=(8, 4))
    
    # plt.pcolormesh(
    #     psi.x,
    #     psi.y,
    #     B_magnitude.T,
    #     shading="auto",
    #     cmap="viridis",
    # )
    
    # plt.colorbar(label=r"$|\mathbf{B}|$ [T]")
    
    # plt.quiver(
    #     X[::skip, ::skip],
    #     Y[::skip, ::skip],
    #     Bx[::skip, ::skip],
    #     By[::skip, ::skip],
    #     color="white",
    #     angles="xy",
    #     scale_units="xy",
    #     scale=None,
    # )
    
    # plt.xlabel("x [m]")
    # plt.ylabel("y [m]")
    # plt.title("Magnetic field magnitude and direction")
    # plt.gca().set_aspect("equal")
    # plt.show()
    
    # plt.figure()
    # plt.imshow(psi.potential,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    # plt.show()
    
    
    

    phi = ElectricPotential(
    Nx=65,
    Ny=65,
    Lx_mm=40.0,
    Ly_mm=60.0,
    E0=1e1,                    # V/m
    gap=0.010,                 # m
    electrode_x_start_mm=-15,   # m
    electrode_x_end_mm= 25,      # m
    )
    
    phi.solve_potential()
    
    ElectricField(phi)
    ElectricField.solve_field()
    ElectricField.show_field()
        



    
    
    psi = MagneticPotential(
        Nx=65,
        Ny=65,
        B0_T=1.0,
        # gap_mm=12.5,
        Lx_mm=12.5,
        Ly_mm=150.0,
        pole_y_start_mm=-75.0,
        pole_y_end_mm=75.0,
    )
    
    psi.solve_potential()
    
    mag = MagneticField(psi)
    mag.solve_field()
    mag.show_field()
    

