# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 10:16:47 2026

@author: TamirCo
"""

def convergence_test (integrator,integrators_name, experiment, engs_MeV,N_steps_vec):
    params = experiment.params
    
    
    
    # colors = {"Boris": "magenta", "Boris_Coll": "red","RK2":"blue", "RK4_Lin": "green", "RK4_Herm": "darkolivegreen", "Euler": "cyan"}
    Interpulations_name = ["Linear","Spline"]
    colors = {"Linear": "darkolivegreen", "Spline":"Olive"}
    R_ends_mm = []
    for i in range (len(Interpulations_name)):
        params["grid_interpulation"] = Interpulations_name[i]
        params = experiment._update_params(params)
        # integrator.mag.interpolation
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        Interpulation_name = Interpulations_name[i]
        R_vec_mm = [];  v_vec_m0s = []; 
        gamma_vec = []; gamma_retio_vec = []; gamma_end_vec = [];
        
        for n in range(len(N_steps_vec)):
            # print (N_steps_vec[n])
            params["N_steps"] = N_steps_vec[n]
            params = experiment._update_params(params)
            integrator._update_params(params)
            R_engs_mm, v_engs_m0s, gamma_engs_vec = basic_simulator (integrator,"RK4_Herm", experiment, engs_MeV,show_tragectories=False)
    
            R_vec_mm.append(R_engs_mm)
            v_vec_m0s.append(v_engs_m0s)
            gamma_vec.append(gamma_engs_vec)
        
        R_end_mm =  [r[-1][-1] for r in R_vec_mm]
        R_ends_mm.append(R_end_mm)
    
    R_ends_mm = np.array(R_ends_mm)
    
    
    
   
    
    
    fig = plt.figure(3)    
    manager = plt.get_current_fig_manager()
    try:
        window = manager.window
     
        window.geometry("1540x900+-10+0")   # WxH+x_offset+y_offset
    except Exception as e:
        print("Could not resize window:", e)

    ax = fig.add_subplot(111)
    def _power_law(x, a, b):
        return a * np.power(x, b)
    
    def _linera_rig(x,a,b):
        return a * x + b
    Y_analitic = integrator._analitic_sol_vel2dist()
    n_fine = np.linspace(N_steps_vec[0],N_steps_vec[-2],1001)
    # power_map = {"Boris": -1,"Boris_Coll": -2 ,"RK2": -2, "RK4_Lin" : -2 ,"RK4_Herm": -4, "Euler": -1}
    power_map = {"Linear": -2, "Spline":-4}
    count = 0
    for n in range (len(Interpulations_name)):
        # params["grid_interpulation"] = Interpulations_name[n]
        # params = experiment._update_params(params)
        # # integrator.mag.interpolation
        # integrator._update_params(params)
        # integrator._update_solution(integrator.solution)
        count += 0.05
        Interpulation_name = Interpulations_name[n]
        Y_min = R_ends_mm[n][-1][1]
        epsilon = np.abs([Y_min - r[1] for r in R_ends_mm[n]])
        # epsilon = np.abs([Y_analitic - r[1] for r in R_ends_mm[n]])
        # print(f"{integrators_name[n]}: epsilon = {epsilon}")
        power = power_map[Interpulation_name]
        log_eps = np.log10(epsilon[:-1])
        log_N = np.log10(N_steps_vec[:-1])
        popt, pcov = curve_fit(_linera_rig, log_N[:], log_eps[:], p0=[power,0],bounds=([-6, -np.inf],[  0,  np.inf]),maxfev=10000)
        a_fit, b_fit = popt
        print(a_fit)
        # eps_fine = _power_law(n_fine,a_fit,b_fit)
        eps_fine = 10**(_linera_rig(np.log10(n_fine[:]),a_fit,b_fit))
        
        color = colors[Interpulation_name]
        N_mid_log = 10**((0.5+count)*(np.log10(N_steps_vec[0]) + np.log10(N_steps_vec[-2])))
        eps_mid_on_fit = 10**(_linera_rig(np.log10(N_mid_log), a_fit, b_fit))
        
        
        ax.text(N_mid_log/2 , eps_mid_on_fit*5,
                rf"$p = {a_fit:.2f}$",
                color=color, fontsize=25, ha='center', va='bottom',
                fontweight='bold')
        label = Interpulation_name
        ax.scatter(N_steps_vec[:-1], epsilon[:-1], color=color, alpha=0.6, s=300,
                   edgecolor='black', linewidth=0.5, label=label)
        ax.plot(n_fine,eps_fine, color=color, linewidth=3)
        
        ax.set_xlabel(r'$N_{steps}$',fontsize=25)
        ax.set_ylabel(r'$\epsilon$',fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.tick_params(axis='both', which='minor', labelsize=18)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(label = "RK4-hermit for Linear Vs Spline grid interpulation",fontsize = 25)
        ax.legend(title='Interpulation',title_fontsize=25,fontsize=18)