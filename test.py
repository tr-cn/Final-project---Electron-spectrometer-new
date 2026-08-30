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

    

def restrict_residual (r_fine):
    # R = 
    return
    


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
        
        residuals = laplace_residuals(phi)
        # eps = err_fun(phi,phi_prev)
        if count % 5 == 0:
            eps = np.max(np.abs(residuals))
    plt.imshow(phi,extent=[-1/2, 1/2, -1/2, 1/2]); plt.colorbar(); 
    plt.title (f'Dericle, Num of iterations:{count} ') 
    plt.show()
    print (count)
    # print(phinew)
        

def main(N = 2**6+1):
    phi0 = np.random.rand(N,N) # [raws,colmns]
    phi0 = apply_boundries (phi0)
    rho = 1 #normlized to epsilon_0
    BC = 0 # dericle
    laplace_func (phi=phi0)
    # BC = 1 # Periodic
    # laplace_func (phi=phi0,BC=BC)

if __name__ == "__main__":
    main()