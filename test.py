import matplotlib.pyplot as plt
import numpy as np
import time

def apply_boundries (phi):
    phi[0,:] = 1; phi[-1,:] = -1; phi[1:-1,0] = 0; phi[1:-1,-1] = 0;
    
    return  phi


def err_fun(phiprev,phinew):
    eps_max = np.max(np.max(np.abs(phinew - phiprev)))\
        
    return eps_max


def gauss_seidel (phi, Nx, Ny): 
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
            R[I,J] = 1/4 * (r [i,j] + r [i,j-1] + r [i-1,j] + r [i-1,-j])
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
    e_fine[0, :] = 0.0
    e_fine[-1, :] = 0.0
    e_fine[:, 0] = 0.0
    e_fine[:, -1] = 0.0

    return e_fine
    

# def restrict_func(phi):
#     Nx, Ny = np.array(phi.shape)
#     ie_swips = 10
#     fe_swips = 100
#     r_finer = laplace_residuals(phi)
#     r_coarser_new = np.copy(r_finer)
#     while Nx>3 or Ny>3:
#         r_finer = np.copy(r_coarser_new)
#         r_coarser = restrict_residual(r_finer)
#         e_coarser = gauss_seidel_residuals(r_coarser,ie_swips)
#         r_coarser_new = r_coarser - e_coarser
#         Nx, Ny = np.array(r_coarser_new.shape)
    
#     # here the grid size is 3
#     e_coarsest = gauss_seidel_residuals(r_coarser,fe_swips);
#     r_coarsest = r_coarser - e_coarsest
 
#     return r_coarsest, e_coarsest


def restrict_func(phi):
    Nx, Ny = np.array(phi.shape)
    ie_swips = 10
    fe_swips = 100
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
    e_coarsest = gauss_seidel_residuals(r_coarser,fe_swips);
    r_coarsest = r_coarser - e_coarsest
    
    r[-1]=(r_coarsest)
    e[-1]=(e_coarsest)
    return r, e





def prolongate_func (e, phi):
    max_Nx, max_Ny = np.array(phi.shape)
    e_coarsest = e[-1]
    Nx, Ny = np.array(e_coarsest.shape)
    count = 0
    while Nx<max_Nx/2 or Ny<max_Ny/2:
        e_finer = prolongate_correction(e[-1-count])    
        e[-1- (count+1)] = e[-1- (count+1)] + e_finer
        count = count + 1
        Nx,Ny = e[-1-count].shape
        
    e_finer = prolongate_correction(e[0])
    phi = phi + e_finer
    
    return phi
        
    




def  laplace_func (phi=[]):
    Nx, Ny = np.array(phi.shape)
    eps = 1;
    treshhold = 1e-6
    count = 0 
    max_round = 10000
    gs_swips = 5 # gauss-seidle swips befor start mutigrid

    
    while eps> treshhold and count<=max_round:
        phi_prev = np.copy(phi)
        count+=1
        for i in range(gs_swips):
            phi = gauss_seidel (phi, Nx, Ny)
        
        
        r, e  = restrict_func(phi)
        phi_prev = np.copy(phi)
        phi = prolongate_func (e, phi)
            
            
        # residuals = laplace_residuals(phi)
        # r_coarse = restrict_residual(residuals)
        
        # e_coarse = gauss_seidel_residuals(r_coarse,ie_swips)
        # e_fine = prolongate_correction(e_coarse)
        
        
        # phi = phi + e_fine
        
        
        
        # eps = err_fun(phi,phi_prev)
        if count % 5 == 0:
            # eps = np.max(np.abs(r[0]))
            eps = np.max(np.abs(phi - phi_prev))
            print (eps)
    plt.imshow(phi,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    plt.title (f'Dericle, Num of iterations:{count} ') 
    plt.show()
    print (count)
    # print(phinew)
        

def main(N = 2**7+1):
    phi0 = np.random.rand(N,N) # [raws,colmns]
    phi0 = apply_boundries (phi0)
    rho = 1 #normlized to epsilon_0
    BC = 0 # dericle
    laplace_func (phi=phi0)
    # BC = 1 # Periodic
    # laplace_func (phi=phi0,BC=BC)

if __name__ == "__main__":
    main()