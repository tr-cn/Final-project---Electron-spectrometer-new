import matplotlib.pyplot as plt
import numpy as np
import time

def apply_boundries (phi):
    phi[0,:] = 0;  phi[-1,:] = 0;
    phi[0,2**4: 2**6] = 1; phi[-1,2**4 :2**6] = -1; phi[1:-1,0] = 0; phi[1:-1,-1] = 0;
    return  phi


def gauss_seidel (phi): 
    Nx, Ny = phi.shape
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            phi[i,j] = 1/4 * (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])
    phi = apply_boundries (phi)
    
    return phi


def laplace_residuals (phi):
    r = np.zeros_like(phi)
    r[1:-1, 1:-1] = (phi[2:, 1:-1] + phi[:-2, 1:-1] + phi[1:-1, 2:] + phi[1:-1, :-2] - 4.0 * phi[1:-1, 1:-1])
    
    return r
 

def restrict_residual (r):
    # r represents the residuals of finer grid
    Nx, Ny = (np.array(r.shape)+ 1)  / 2
    Nx = int(Nx); Ny = int(Ny)
    R = np.zeros([Nx,Ny])
    for I in range (1,Nx-1):
        for J in range (1,Ny-1):
            i = 2*I
            j = 2*J
            R[I,J] = 1/4 * (r [i,j] + r [i,j-1] + r [i-1,j] + r [i-1,j-1])
    r[0,:] = 0; r[-1,:] = 0; r[:,0] = 0; r[:,-1] = 0;
    
    return R
    

def gauss_seidel_residuals(R,N):
    
    # r : residuals
    # N : number of iterations
    Nx,Ny = R.shape
    
    e = np.zeros_like(R)
    for n in range (N):
        for i in range (1,Nx-1):
            for j in range (1,Ny-1):
                e[i,j] = 1/4 * (e[i+1,j] + e[i-1,j] + e[i,j+1] + e[i,j-1] + R[i,j])
        
    e[0, :] = 0.0; e[-1, :] = 0.0; e[:, 0] = 0.0; e[:, -1] = 0.0;
    
    return e
        
    
def prolongate_correction(e_coarse):
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
    




def restrict_func(phi):
    Nx, Ny = phi.shape
    ie_swips = 4
    fe_swips = 10
    r = []; e = []
        
    r.append(laplace_residuals(phi))
    while Nx>3 or Ny>3:
        r_finer = r[-1]
        r_coarser = restrict_residual(r_finer)
        e_coarser = gauss_seidel_residuals(r_coarser,ie_swips)
        r_coarser_new = r_coarser - laplace_residuals(e_coarser)
       
        r.append(r_coarser_new)
        e.append(e_coarser)
        Nx, Ny = np.array(r_coarser_new.shape)
    
    # here the grid size is 3
    e_coarsest = gauss_seidel_residuals(r[-1],fe_swips);
    r_coarsest = r_coarser - laplace_residuals(e_coarsest)
    
    r[-1]=(r_coarsest)
    e[-1]=(e_coarsest)
    return r, e

def phi_smother(phi,gs_swips):
    
    for i in range(gs_swips):
        phi = gauss_seidel(phi)
    phi = apply_boundries(phi)
    return phi
    

def prolongate_func(e, phi):
    # Go from the coarsest correction upward
    for level in range(len(e) - 2, -1, -1):
        finer_correction = prolongate_correction(e[level + 1])

        # Add the correction to the current level
        e[level] = e[level] + finer_correction

    # Transfer the finest correction to phi
    phi = phi + prolongate_correction(e[0])

    # Restore the physical boundary conditions
    phi = apply_boundries(phi)
   

    return phi




def  laplace_func (phi):
    Nx, Ny = np.array(phi.shape)
    eps = 1;
    treshhold = 1e-6
    count = 0 
    max_round = 10000
    gs_swips = 5 # gauss-seidle swips befor start mutigrid

    
    while eps > treshhold and count<=max_round:
        phi_prev = np.copy(phi)
        count+=1
        
        phi_smother(phi,gs_swips)
        r, e  = restrict_func(phi)
        phi = prolongate_func (e, phi)
        phi = phi_smother(phi,gs_swips)
      
        if count % 5 == 0:
            eps = np.max(np.abs(laplace_residuals(phi)))
            print (eps)
    plt.imshow(phi,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    plt.title (f'Dericle, Num of iterations:{count} ') 
    plt.show()
    print (count)
    # print(phinew)
    
    return phi
        








def get_elctric_potential(Nx = 2**6+1, Ny = 2**6+1):
    phi0 = np.random.rand(Nx,Ny) # [raws,colmns]
    phi0 = apply_boundries (phi0)
    phi = laplace_func (phi0)
    return phi
    

def get_magnetic_potential(Nx = 2**6+1, Ny = 2**6+1):
    mu0 = 4.0 * np.pi * 1e-7
    B0_T = 0.5       # Tesla
    D_m = 0.01       # m
    delta_psi_A = B0_T * D_m / mu0
    
    psi0 = np.random.rand(Nx,Ny) # [raws,colmns]
    psi0 = apply_boundries (psi0)
    psi0 = psi0 * delta_psi_A
    psi0 = laplace_func (psi0)
    
    
    return


# def main():
#     Nx = 2**6+1; Ny = 2**6+1;
#     phi = get_elctric_potential(Nx =Nx, Ny = Ny)
       
    
#     return phi
    
    # BC = 1 # Periodic
    # laplace_func (phi=phi0,BC=BC)

if __name__ == "__main__":
    plt.close('all')
    # phi = main()
    Nx = 2**6+1; Ny = 2**6+1;
    # phi = get_elctric_potential(Nx = Nx, Ny = Ny)
    psi = get_magnetic_potential(Nx = Nx, Ny = Ny)