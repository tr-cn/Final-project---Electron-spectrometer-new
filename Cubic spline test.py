import matplotlib.pyplot as plt
import numpy as np
import time
plt.close("all")

global Delta_x, S_intervals

def func(x,lamb):
    return np.sin((x-x**2))/x + lamb*np.random.uniform(low=-1.0, high=1.0, size=x.shape)

def P_tag_tag_fun (vec,f): # Sloving the A p'' = b set of linear eqautions
    global Delta_x, interval
    
    interval = len(vec)-1 # number of interveal
    Delta_x = (vec[-1] - vec[0])/interval # assuming eqaul intervals
    M_size = interval - 1 # Matrix size
    
    # Main diagonal: [4*dx, 4*dx, ...]
    main_diag = 4 * Delta_x * np.ones(M_size)

    # Off-diagonals (sub and super): [dx, dx, ...]
    off_diag = Delta_x * np.ones(M_size - 1)

    # Build the tridiagonal matrix A
    A = np.diag(main_diag, k=0) +\
        np.diag(off_diag, k=1) + \
        np.diag(off_diag, k=-1)



    # Calculate b vector
    b_inner = f[:-2] - 2 * f[1:-1] + f[2:] 
    b = (6 / Delta_x) * b_inner

    # find p''
    p_tag_tag = np.linalg.solve(A, b)
    p_tag_tag = np.append(np.array(0), p_tag_tag)
    p_tag_tag = np.append(p_tag_tag, np.array(0)) 

    return (p_tag_tag)

def Params_finder (p_tag_tag): # Finding the alph, beta, gammae, eta coefficents
    global Delta_x
    alpha = np.array([]); beta = alpha; gamma = beta; etha = gamma;
    
    for i in range(len(p_tag_tag)-1):
    #print(alpha.shape)
        alpha = np.append(alpha,(p_tag_tag[i+1]/(6*Delta_x)))
        beta = np.append(beta,(-p_tag_tag[i]/(6*Delta_x)))
        gamma = np.append(gamma,(-p_tag_tag[i+1]*Delta_x**2+6*f[i+1])/(6*Delta_x))
        etha = np.append(etha,(p_tag_tag[i]*Delta_x**2-6*f[i])/(6*Delta_x) )

    return alpha, beta, gamma, etha
    

    
    
def f_spline(vec, alpha, beta, gamma, etha): # Genereating the splined function
    global S_intervals
    for s in range(len(vec)-1):
    
        x = np.linspace(vec[s],vec[s+1],S_intervals) # making each interval been "continues"
        
        a = np.array( alpha[s] * (x-vec[s])**3 )
        b = np.array( beta[s]  * (x-vec[s+1])**3 )
        g = np.array( gamma[s] * (x-vec[s]) )
        e = np.array( etha[s]  * (x-vec[s+1]) )
        v_tot = a+b+g+e
        plt.plot(x,v_tot)
        plt.show()
        
    
    
def f_tag_spline (vec, alpha, beta, gamma, etha): # Finding the first derivative of the spline
    global S_intervals
    fig = plt.figure
    for s in range(len(vec)-1):
        
        x = np.linspace(vec[s],vec[s+1],S_intervals)
        
        a_tag = np.array( alpha[s] * 3* (x-vec[s])**2 )
        b_tag = np.array( beta[s]  * 3*(x-vec[s+1])**2 )
        g_tag = np.array( gamma[s] )
        e_tag = np.array( etha[s] )    
        
        v_tot_tag = a_tag+b_tag+g_tag+e_tag
        plt.plot(x,v_tot_tag)
        plt.show()
        

S_intervals = 1001
vec = np.linspace(0.5,10,500)
lambda_vec = np.array([0.1,0.05,0.01,0.005,0])#
for l in lambda_vec:
    f = func(vec,l)  


    p_tag_tag = P_tag_tag_fun (vec,f)
    #print(len(P_tag_tag))
    alpha, beta, gamma, etha = Params_finder (p_tag_tag)
    #print(len(alpha))
    f_spline(vec, alpha, beta, gamma, etha); plt.plot(vec,f,color = 'black')
    plt.title(f' f vs f_spline: lambda = {l}'); plt.show()
#     print(len(f_tag))
    f_tag = (f[2:-1] - f[0:-3])/(2*Delta_x)
    f_tag_spline (vec[0:-1], alpha, beta, gamma, etha); plt.plot(vec[1:-2],f_tag,color = 'black')
    plt.title( f' f_tag vs f_tag_spline:lambda = {l}');
    plt.show()
    
    # %%
    
    def _spine_2D(self, R_m):
        def __P_tag_tag_fun (vec,f): # Sloving the A p'' = b set of linear eqautions
            
            
            interval = len(vec)-1 # number of interveal
            Delta_x = (vec[-1] - vec[0])/interval # assuming eqaul intervals
            M_size = interval - 1 # Matrix size
            
            # Main diagonal: [4*dx, 4*dx, ...]
            main_diag = 4 * Delta_x * np.ones(M_size)

            # Off-diagonals (sub and super): [dx, dx, ...]
            off_diag = Delta_x * np.ones(M_size - 1)

            # Build the tridiagonal matrix A
            A = np.diag(main_diag, k=0) +\
                np.diag(off_diag, k=1) + \
                np.diag(off_diag, k=-1)



            # Calculate b vector
            b_inner = f[:-2] - 2 * f[1:-1] + f[2:] 
            b = (6 / Delta_x) * b_inner

            # find p''
            p_tag_tag = np.linalg.solve(A, b)
            p_tag_tag = np.append(np.array(0), p_tag_tag)
            p_tag_tag = np.append(p_tag_tag, np.array(0)) 

            return p_tag_tag, interval
        
        def __f_tag_spline (vec, alpha, beta, gamma, etha,interval): # Finding the first derivative of the spline
            
            for s in range(len(vec)-1):
            
                x = np.linspace(vec[s],vec[s+1],S_intervals)
                
                a_tag = np.array( alpha[s] * 3* (x-vec[s])**2 )
                b_tag = np.array( beta[s]  * 3*(x-vec[s+1])**2 )
                g_tag = np.array( gamma[s] )
                e_tag = np.array( etha[s] )    
                
                v_tot_tag = a_tag+b_tag+g_tag+e_tag
                
        def __f_spline(vec, alpha, beta, gamma, etha,interval): # Genereating the splined function
            
            for s in range(len(vec)-1):
            
                x = np.linspace(vec[s],vec[s+1],S_intervals) # making each interval been "continues"
                
                a = np.array( alpha[s] * (x-vec[s])**3 )
                b = np.array( beta[s]  * (x-vec[s+1])**3 )
                g = np.array( gamma[s] * (x-vec[s]) )
                e = np.array( etha[s]  * (x-vec[s+1]) )
                v_tot = a+b+g+e
                
                
        
        p_tag_tag = P_tag_tag_fun (vec,f)
        alpha, beta, gamma, etha = Params_finder (p_tag_tag)
        f_spline(vec, alpha, beta, gamma, etha); 

    #     print(len(f_tag))
        f_tag = (f[2:-1] - f[0:-3])/(2*Delta_x)
        f_tag_spline (vec[0:-1], alpha, beta, gamma, etha)
        
        
    
    
    
