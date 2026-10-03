#from IPython import get_ipython
#import matplotlib.pyplot as plt
#plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');


# %%
import matplotlib.pyplot as plt

import numpy as np
from types import SimpleNamespace

plt.close('all')

from spectrometer.geometry import Spectrometer_Budy
from spectrometer.Experiment import Experiment
from spectrometer.integrators import Integrators
# import spectrometer.integrators as integ
import spectrometer.plot_result as pr
import spectrometer.analitic_fields as analitic_fields
from spectrometer.simulators import *
import copy




if __name__ == "__main__":
   
    
    me_kg =9.109*1e-31 #9.109*1e-31
    R0_mm  = np.array([0,-10,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = 1e8

    yoke = 0
    shield_mm = 5
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.05
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =    "Numeric Field"  #  "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**4 + 1;
    Ny_p = 2**4 + 1;
    
    grid_interpulation = "Linear" # "Spline"#
    # spec = Spectrometer_Budy(shield_mm = shield_mm, yoke=yoke)
    # fig,ax = spec._draw_spec()
    
    experiment = Experiment(
                            height_mm = height_mm, width_mm = width_mm, depth_mm = depth_mm, shield_mm = shield_mm, pinhole_dia_mm = pinhole_dia_mm, yoke = yoke, # Spectrometere and simulation border
                            R0_mm= R0_mm, q_eng_MeV = e_eng_MeV, m_kg = me_kg, q_C = e_C, # Particle parameters
                            Bx0_T = Bx0_T, Ex0_Vm = Ex0_Vm, fringe = fringe, sharp_edge = sharp_edge, # Fildes parameters
                            CFL = CFL, N_steps = N_steps, solution = solution, # Simulation resolutions
                            Nx_p = Nx_p, Ny_p = Ny_p, grid_interpulation = grid_interpulation, # Potential solving
                            )
       
    experiment._evaluate_exp_paramas()
    params = experiment.params
   
    # %%    #_________________Integration environment intitiation__________________
    if True:#False:
        params["solution"] = "Analitic field"
    
        integrator = Integrators(params)
    
    if False:#True:
        params["solution"] = "Numeric Field"
        integrator = Integrators(params)
        
        # integrator.psi._show_potential()
        # integrator.mag._show_field()
        # integrator.phi._show_potential()
        # integrator.Ele._show_field()
        
    params_beckup = copy.deepcopy(params)
    #%%    Fully analitic assuming no velocity in Z direction, no fringe fields
    if True:
        params = copy.deepcopy(params_beckup)
        params["Ex0_Vm"] = 5e7
        params["R0_mm"] = np.array([0,-4,0])
        
        params = experiment._update_params(params)
        # integrator.mag.interpolation
        params = integrator._update_params(params)
        e_engs_vec_MeV = np.array([ [0,i,0] for i in np.linspace(1,15,1001)])
        m_kg_vec  =[9.109*1e-31,9.1e-30]# [9.1e-30]#[9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,-1.602*1e-18] #[-1.602*1e-18]#[-1.602*1e-19,-1.602*1e-18]
        # params["q_eng_MeV"] = np.array([10,10,0])
        Fully_analitic(integrator,experiment,e_engs_vec_MeV,m_kg_vec,q_vec_C)
        
        
        
        
        
    # %%    #___________________________Basic simuations_______________________
    
    
    
    if False:
        
        e_engs_MeV = np.array([ [0,i,0] for i in np.linspace(1,15,15)])
        # e_engs_MeV = np.array([-0.005,2,0.005])# When using this, it is better to work with velocity in the y driection so the energy will not exceed the speed of light
        
        integrator_name = "Boris"#"Euler"#"RK4"#"RK2"
        basic_simulator(integrator,integrator_name, experiment, e_engs_MeV,show_tragectories = True)
    

    # %%    #___________________________Converssion Tests______________________
    if False:
        params = copy.deepcopy(params_beckup)
        integrators_name =["RK4_Herm"]#, "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]#["RK4_Herm"]#
        params["fringe"] = 1
        params["Ex0_Vm"] = Ex0_Vm
        params["R0_mm"] = np.array([0,-4,0])
        params["grid_interpulation"] = "Spline"#"Spline"# "Linear"
        params = experiment._update_params(params)
        # integrator.mag.interpolation
        params = integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        e_engs_MeV = np.array([0, 3, 0])
        # N_steps_vec = [1e2,5e2,1e3,1e4]
# %%
        N_steps_vec = [25,50,100,200,400,800,1600,3200,6400]

        convergence_test (integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
          
# %%    # ___________________________Boris Tests_______________________________  
    if False:
        def _compare_boris_vs_rk4(integrator, experiment, N_steps_test=5000):
            """
            מריץ Boris ו-RK4 עם אותם תנאים, ומדפיס את המסלול השלם של שניהם
            כדי להשוות איפה הם נבדלים.
            """
            params = experiment.params
            params["N_steps"] = N_steps_test
            params = experiment._update_params(params)
            integrator._update_params(params)
        
            R_mm_rk4, v_rk4, gamma_rk4 = integrator._RK4()
            R_vec_rk4 = [r.copy() for r in integrator.R_vec_m]
        
            integrator._update_params(params)
            R_mm_boris, v_boris, gamma_boris = integrator._Boris_pusher()
            R_vec_boris = [r.copy() for r in integrator.R_vec_m]
        
            print(f"RK4 steps: {len(R_vec_rk4)}, Boris steps: {len(R_vec_boris)}")
            print(f"\nRK4 final R (m):   {R_vec_rk4[-1]}")
            print(f"Boris final R (m): {R_vec_boris[-1]}")
            print(f"Difference: {R_vec_rk4[-1] - R_vec_boris[-1]}")
        
            n_compare = min(len(R_vec_rk4), len(R_vec_boris), 10)
            print(f"\nFirst {n_compare} steps comparison (y-coordinate):")
            for i in range(n_compare):
                print(f"  step {i}: RK4_y={R_vec_rk4[i][1]:.8e}, Boris_y={R_vec_boris[i][1]:.8e}, "
                      f"diff={abs(R_vec_rk4[i][1]-R_vec_boris[i][1]):.3e}")
        
            print(f"\nLast {n_compare} steps comparison (y-coordinate):")
            for i in range(-n_compare, 0):
                print(f"  step {i}: RK4_y={R_vec_rk4[i][1]:.8e}, Boris_y={R_vec_boris[i][1]:.8e}, "
                      f"diff={abs(R_vec_rk4[i][1]-R_vec_boris[i][1]):.3e}")
        
            return R_vec_rk4, R_vec_boris
        def _boris_convergence_order(integrator, experiment, N_list=None):
            """
                                מריץ Boris בכמה ערכי N_steps ומחשב את סדר-ההתכנסות ש            לפי Y_exit(N) מול Y_exit(2N) -- Richardson extrapolation.
            """
            if N_list is None:
                N_list = [1000, 2000, 4000, 8000, 16000]
            
            Y_exit_list = []
            for N in N_list:
                params = experiment.params
                params["N_steps"] = N
                params = experiment._update_params(params)
                integrator._update_params(params)
            
                R_mm, v_m0s, gamma_vec = integrator._Boris_pusher()
                Y_exit_mm = R_mm[-1][1]
                Y_exit_list.append(Y_exit_mm)
                print(f"N_steps={N:>7d}  Y_exit={Y_exit_mm:.10f} mm  actual_steps={len(R_mm)}")
            
            Y_exit_arr = np.array(Y_exit_list)
            
            print("\n--- Richardson extrapolation (epsilon = |Y(N) - Y(2N)|) ---")
            eps_list = []
            for i in range(len(N_list)-1):
                eps = abs(Y_exit_arr[i] - Y_exit_arr[i+1])
                eps_list.append(eps)
                print(f"N={N_list[i]:>7d} vs N={N_list[i+1]:>7d}:  epsilon = {eps:.6e}")
            
            print("\n--- Observed convergence order p (eps_i / eps_{i+1} = 2^p) ---")
            orders = []
            for i in range(len(eps_list)-1):
                if eps_list[i+1] == 0:
                    continue
                ratio = eps_list[i] / eps_list[i+1]
                p = np.log2(ratio)
                orders.append(p)
                print(f"eps({N_list[i]}->{N_list[i+1]}) / eps({N_list[i+1]}->{N_list[i+2]}) = {ratio:.4f}  =>  order p = {p:.4f}")
            
            if orders:
                print(f"\nAverage observed order: {np.mean(orders):.4f}")
            
            return N_list, Y_exit_list, eps_list, orders
        def boris_order_no_exit(integrator, experiment, N_list=None, n_steps_fixed=200):
            """
                    בודק את סדר-ההתכנסות של Boris בלי מנגנון ה-exit-crossing בכלל        מריץ מספר-צעדים *קבוע* (לא עד שחוצה גבול), ובודק אם Y מתכנס בסדר-1 או סד        כש-N_steps (צפיפות הצעדים הפנימית) משתנה.
            """
            import numpy as np
            if N_list is None:
                N_list = [1000, 2000, 4000, 8000, 16000, 32000]
            
            Y_list = []
            for N in N_list:
                params = experiment.params
                params["N_steps"] = N
                params = experiment._update_params(params)
                integrator._update_params(params)
            
                dt_s = integrator.dt_s
                R_current_m = integrator.R0_m
                E_Vm = integrator._get_electric_field(R_current_m)
                from spectrometer.physics import get_electric_acceleration
                dv_dt_0 = get_electric_acceleration(integrator.q_C, integrator.gamma0, integrator.m_kg, E_Vm, integrator.v0_m0s, B_T=None)
                v_current_m0s = integrator.v0_m0s - 0.5*dv_dt_0*dt_s
                gamma_current = integrator.gamma0
            
                from spectrometer.physics import get_magnetic_rotation
                from spectrometer.integrators import vel2gamma
            
                for _ in range(n_steps_fixed):
                    B_T = integrator._get_magnetic_field(R_current_m)
                    E_Vm_i = integrator._get_electric_field(R_current_m)
            
                    dv1 = get_electric_acceleration(integrator.q_C, gamma_current, integrator.m_kg, E_Vm_i, v_current_m0s, B_T=None) * dt_s/2
                    v1 = v_current_m0s + dv1
                    gamma_m1 = vel2gamma(v1)
            
                    v2 = get_magnetic_rotation(integrator.q_C, gamma_m1, integrator.m_kg, E_Vm_i, v1, B_T, dt_s)
                    gamma_m2 = vel2gamma(v2)
            
                    v_current_m0s = v2 + get_electric_acceleration(integrator.q_C, gamma_m2, integrator.m_kg, E_Vm_i, v2, B_T=None)*dt_s/2
                    R_current_m = R_current_m + v_current_m0s*dt_s
                    gamma_current = vel2gamma(v_current_m0s)
            
                Y_fixed_time_mm = R_current_m[1]*1e3
                physical_time_s = n_steps_fixed * dt_s
                Y_list.append(Y_fixed_time_mm)
                print(f"N_steps={N:>7d}  dt_s={dt_s:.6e}  physical_t={physical_time_s:.6e}s  Y(after {n_steps_fixed} steps)={Y_fixed_time_mm:.10f} mm")
            
            print("\n--- epsilon and order, FIXED physical time, no exit-crossing ---")
            import numpy as np
            Y_arr = np.array(Y_list)
            eps_list = [abs(Y_arr[i]-Y_arr[i+1]) for i in range(len(Y_arr)-1)]
            for i,eps in enumerate(eps_list):
                print(f"N={N_list[i]}->{N_list[i+1]}: epsilon={eps:.6e}")
            orders=[]
            for i in range(len(eps_list)-1):
                if eps_list[i+1]==0: continue
                p = np.log2(eps_list[i]/eps_list[i+1])
                orders.append(p)
                print(f"order p = {p:.4f}")
            if orders:
                print(f"\nAverage order (no exit-crossing): {np.mean(orders):.4f}")
            return N_list, Y_list, eps_list, orders   
        def boris_order_fixed_physical_time(integrator, experiment, N_list=None, t_fixed_s=None):
            """
            בודק את סדר-ההתכנסות של Boris בלי מנגנון ה-exit-crossing,
            עד *אותו זמן-פיזיקלי קבוע* בכל הריצות (לא מספר-צעדים קבוע!).
            """
            import numpy as np
            from spectrometer.physics import get_electric_acceleration, get_magnetic_rotation
            from spectrometer.integrators import vel2gamma
        
            if N_list is None:
                N_list = [1000, 2000, 4000, 8000, 16000, 32000]
        
            if t_fixed_s is None:
                params = experiment.params
                params["N_steps"] = N_list[0]
                params = experiment._update_params(params)
                integrator._update_params(params)
                t_fixed_s = 200 * integrator.dt_s
                print(f"t_fixed_s = {t_fixed_s:.6e} s (קבוע עבור כל הריצות)")
        
            Y_list = []
            for N in N_list:
                params = experiment.params
                params["N_steps"] = N
                params = experiment._update_params(params)
                integrator._update_params(params)
        
                dt_s = integrator.dt_s
                n_steps_for_this_N = int(round(t_fixed_s / dt_s))
        
                R_current_m = integrator.R0_m
                E_Vm = integrator._get_electric_field(R_current_m)
                dv_dt_0 = get_electric_acceleration(integrator.q_C, integrator.gamma0, integrator.m_kg, E_Vm, integrator.v0_m0s, B_T=None)
                v_current_m0s = integrator.v0_m0s - 0.5*dv_dt_0*dt_s
                gamma_current = integrator.gamma0
        
                for _ in range(n_steps_for_this_N):
                    B_T = integrator._get_magnetic_field(R_current_m)
                    E_Vm_i = integrator._get_electric_field(R_current_m)
        
                    dv1 = get_electric_acceleration(integrator.q_C, gamma_current, integrator.m_kg, E_Vm_i, v_current_m0s, B_T=None) * dt_s/2
                    v1 = v_current_m0s + dv1
                    gamma_m1 = vel2gamma(v1)
        
                    v2 = get_magnetic_rotation(integrator.q_C, gamma_m1, integrator.m_kg, E_Vm_i, v1, B_T, dt_s)
                    gamma_m2 = vel2gamma(v2)
        
                    v_current_m0s = v2 + get_electric_acceleration(integrator.q_C, gamma_m2, integrator.m_kg, E_Vm_i, v2, B_T=None)*dt_s/2
                    R_current_m = R_current_m + v_current_m0s*dt_s
                    gamma_current = vel2gamma(v_current_m0s)
        
                Y_mm = R_current_m[1]*1e3
                Y_list.append(Y_mm)
                actual_t = n_steps_for_this_N * dt_s
                print(f"N_steps={N:>7d}  dt_s={dt_s:.6e}  n_inner_steps={n_steps_for_this_N:>6d}  actual_t={actual_t:.6e}s  Y={Y_mm:.10f} mm")
        
            Y_arr = np.array(Y_list)
            print("\n--- epsilon and order, FIXED physical time, no exit-crossing ---")
            eps_list = [abs(Y_arr[i]-Y_arr[i+1]) for i in range(len(Y_arr)-1)]
            for i,eps in enumerate(eps_list):
                print(f"N={N_list[i]}->{N_list[i+1]}: epsilon={eps:.6e}")
            orders=[]
            for i in range(len(eps_list)-1):
                if eps_list[i+1]==0: continue
                p = np.log2(eps_list[i]/eps_list[i+1])
                orders.append(p)
                print(f"order p = {p:.4f}")
            if orders:
                print(f"\nAverage order (fixed physical time, no exit-crossing): {np.mean(orders):.4f}")
            return N_list, Y_list, eps_list, orders
        def boris_order_fixed_time_correct_init(integrator, experiment, N_list=None, t_fixed_s=None):
            """
            זהה לקודם, אבל מתקן את אתחול v_minus_half כך שיכלול rotation מגנטי אמיתי
            (half-step-back מלא, לא רק E-kick) -- לבדוק אם זה מתקן את הסדר ל-2.
            """
            import numpy as np
            from spectrometer.physics import get_electric_acceleration, get_magnetic_rotation
            from spectrometer.integrators import vel2gamma
        
            if N_list is None:
                N_list = [1000, 2000, 4000, 8000, 16000, 32000]
        
            if t_fixed_s is None:
                params = experiment.params
                params["N_steps"] = N_list[0]
                params = experiment._update_params(params)
                integrator._update_params(params)
                t_fixed_s = 200 * integrator.dt_s
        
            Y_list = []
            for N in N_list:
                params = experiment.params
                params["N_steps"] = N
                params = experiment._update_params(params)
                integrator._update_params(params)
        
                dt_s = integrator.dt_s
                n_steps_for_this_N = int(round(t_fixed_s / dt_s))
        
                R_current_m = integrator.R0_m
                B_T0 = integrator._get_magnetic_field(R_current_m)
                E_Vm0 = integrator._get_electric_field(R_current_m)
        
                # half-step-back תקין: E-kick(-dt/4) + B-rotation(-dt/2) + E-kick(-dt/4) בקירוב,
                # או לפחות לכלול סיבוב-B חלקי לאחור, לא רק E
                dv_E0 = get_electric_acceleration(integrator.q_C, integrator.gamma0, integrator.m_kg, E_Vm0, integrator.v0_m0s, B_T=None)
                v_temp = integrator.v0_m0s - 0.5*dv_E0*dt_s
                gamma_temp = vel2gamma(v_temp)
                # סיבוב-B לאחור (dt שלילי, half-step-back מלא)
                v_minus_half = get_magnetic_rotation(integrator.q_C, gamma_temp, integrator.m_kg, E_Vm0, v_temp, B_T0, -dt_s)
        
                v_current_m0s = v_minus_half
                gamma_current = vel2gamma(v_current_m0s)
        
                for _ in range(n_steps_for_this_N):
                    B_T = integrator._get_magnetic_field(R_current_m)
                    E_Vm_i = integrator._get_electric_field(R_current_m)
        
                    dv1 = get_electric_acceleration(integrator.q_C, gamma_current, integrator.m_kg, E_Vm_i, v_current_m0s, B_T=None) * dt_s/2
                    v1 = v_current_m0s + dv1
                    gamma_m1 = vel2gamma(v1)
        
                    v2 = get_magnetic_rotation(integrator.q_C, gamma_m1, integrator.m_kg, E_Vm_i, v1, B_T, dt_s)
                    gamma_m2 = vel2gamma(v2)
        
                    v_current_m0s = v2 + get_electric_acceleration(integrator.q_C, gamma_m2, integrator.m_kg, E_Vm_i, v2, B_T=None)*dt_s/2
                    R_current_m = R_current_m + v_current_m0s*dt_s
                    gamma_current = vel2gamma(v_current_m0s)
        
                Y_mm = R_current_m[1]*1e3
                Y_list.append(Y_mm)
                print(f"N_steps={N:>7d}  dt_s={dt_s:.6e}  Y={Y_mm:.10f} mm")
        
            Y_arr = np.array(Y_list)
            print("\n--- epsilon and order, corrected init, fixed physical time ---")
            eps_list = [abs(Y_arr[i]-Y_arr[i+1]) for i in range(len(Y_arr)-1)]
            for i,eps in enumerate(eps_list):
                print(f"N={N_list[i]}->{N_list[i+1]}: epsilon={eps:.6e}")
            orders=[]
            for i in range(len(eps_list)-1):
                if eps_list[i+1]==0: continue
                p = np.log2(eps_list[i]/eps_list[i+1])
                orders.append(p)
                print(f"order p = {p:.4f}")
            if orders:
                print(f"\nAverage order (corrected init): {np.mean(orders):.4f}")
            return N_list, Y_list, eps_list, orders
        
      
        # R_vec_rk4, R_vec_boris = _compare_boris_vs_rk4(integrator, experiment, N_steps_test=5000)
        
        # N_list, Y_exit_list, eps_list, orders = _boris_convergence_order(
                # integrator, experiment, N_list=[1000, 2000, 4000, 8000, 16000, 32000]                )
        
        # N_list2, Y_list2, eps_list2, orders2 = boris_order_no_exit(integrator, experiment, N_list=[1000,2000,4000,8000,16000,32000], n_steps_fixed=200)
        
        # N_list3, Y_list3, eps_list3, orders3 = boris_order_fixed_physical_time(integrator, experiment, N_list=[1000,2000,4000,8000,16000,32000])
        
        N_list4, Y_list4, eps_list4, orders4 = boris_order_fixed_time_correct_init(integrator, experiment, N_list=[1000,2000,4000,8000,16000,32000])
            # print(R_vec_rk4)
            # print(R_vec_boris)

    # %%    # _____________Energy Vs N_steps Vs Integrator_____________________
    if False:
        params = copy.deepcopy(params_beckup)
        integrators_name =["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 0
        params["R0_mm"] = np.array([0,0,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        e_engs_MeV = np.array([0,5,0])
        N_steps_vec = [1e2]#,4e3,5e3,1e4,3e4,5e4]
         
        energy_vs_N_steps_for_different_integrators(integrator,integrators_name, experiment, e_engs_MeV,N_steps_vec)
        
        
    # %%# _____________________Non_collimataed_beams___________________________
    
    if False:
        params = copy.deepcopy(params_beckup)
        integrators_name = ["RK4_Herm"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 1
        params["Ex0_Vm"] = 1e-7
        params["R0_mm"] = np.array([0,-3,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        # integrator._update_solution(integrator.solution)
        
        e_engs_MeV_scalar = 4
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Exponential"#"Uniform" #"Gaussian" #"Exponential"
        Non_collimataed_beams(integrator,integrators_name, experiment,N_steps_vec, e_engs_MeV_scalar, number_of_experiments, angular_distribution,angular_scale=0.03, show_tragectories = False, show_map = True)
        
    # %% # ______________Real_beams - spectral distribution__________________ 
    if False:
        
        integrator_name = ["RK4_Lin"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 0
        params["Ex0_Vm"] = 1e6
        params["R0_mm"] = np.array([0,-5,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        e_engs_MeV_scalar = 4
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential" 
        angular_scale = 0.005
        energy_scale = 0.1
        real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
                   N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
    # %% # _______________________TNSA - experiment_________________________
    if True:
        params = copy.deepcopy(params_beckup)
        params = experiment.params 
        integrator_name = ["RK4_Lin"]#["Euler", "Boris","Boris_Coll","RK2","RK4_Lin", "RK4_Herm"]
        params["fringe"] = 1
        params["Ex0_Vm"] = 5e7
        params["R0_mm"] = np.array([0,-4,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        integrator._update_solution(integrator.solution)
        
        
        
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Gaussian"#"Gaussian" #"None"#,"Uniform" #"Gaussian" #"Exponential"
        energy_distribution =  "DLA" #"None"#"Uniform" #"Gaussian" #"Exponential"#"DLA"
        angular_scale = 0.005
        energy_scale = 0.1
        e_engs_MeV_scalars = [3,0.1]
        mq_vec_kg = [9.109*1e-31,9.1e-30]
        q_vec_C = [-1.602*1e-19,-1.602*1e-18]
 
        TNSA (integrator = integrator,integrator_name =integrator_name, experiment =experiment, mq_vec_kg = mq_vec_kg, q_vec_C = q_vec_C,
                   N_steps = N_steps, total_energy_scalars_MeV = e_engs_MeV_scalars, energy_distribution= energy_distribution, energy_scale = energy_scale,
                   number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
                   angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
        
        
            
#%%        
    # if True:
    #     integrator_name = ["Boris"]#["Boris","RK2","RK4","Euler"]
    #     e_engs_MeV_scalar = 4
    #     N_steps_vec = 1e3#,4e3,5e3,1e4,3e4,5e4]
    #     number_of_experiments = 60
    #     angular_distribution = "Uniform" #"Gaussian" #"Exponential"
    #     energy_distribution =  "Uniform" #"Gaussian" #"Exponential"
    #     angular_scale = 0.03
    #     energy_scale = 0.1
    #     real_beams(integrator = integrator,integrator_name =integrator_name, experiment =experiment,
    #                N_steps = N_steps, total_energy_scalar_MeV = e_engs_MeV_scalar, energy_distribution= energy_distribution, energy_scale = energy_scale,
    #                number_of_experiments = number_of_experiments, angular_distribution = angular_distribution,
    #                angular_scale=angular_scale, show_tragectories = False, show_map = True)
        
        
        
        
# def dNdE_shape(E, Te):
#     return np.sqrt(E/np.pi) * Te**(-1.5) * np.exp(-E/Te)

# Te = 3.0   # טמפרטורת-האלקטרונים, MeV
# E_peak = Te / 2                              # שיא הפונקציה (נגזרת אנליטית)
# f_max = dNdE_shape(E_peak, Te)
# E_max_range = 15 * Te                        # תחום גזירה, כולל ה"זנב" האקספוננציאלי

# def rejection_sample_energy(N, Te, E_max_range, f_max):
#     samples = []
#     while len(samples) < N:
#         E_candidates = np.random.uniform(0, E_max_range, N)
#         y_candidates = np.random.uniform(0, f_max, N)
#         accepted = E_candidates[y_candidates <= dNdE_shape(E_candidates, Te)]
#         samples.extend(accepted.tolist())
#     return np.array(samples[:N])

# energies_MeV = rejection_sample_energy(N=1000, Te=Te, E_max_range=E_max_range, f_max=f_max)
    
    
    
    # simulation_engs_collimated(integrator_anal, experiment, e_engs_MeV,ax)
    
    
    
        
    
    
    # integrator.solution = "Numeric Field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._RK4 ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # integrator.solution = "Analitic field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._RK4 ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    
    
    # integrator.solution = "Numeric Field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._Boris_pusher ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # integrator.solution = "Analitic field"
    # R_vec_mm,v_vec_m0s,gamma_vec = integrator._Boris_pusher ()
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    
    
    
    
    
    
    
    # Z_mm = integ.analitic_sol_vel2dist (e_eng_MeV[1], me_kg, e_C, height_mm, B_T[0])
    # print(Z_mm)
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.euler (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.RK2 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.RK4 (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    # R_vec_mm,v_vec_m0s,gamma_vec = integ.Boris_pusher (e_eng_MeV, me_kg, e_C, height_mm, width_mm, depth_mm, B_T, E_V0m, R0_mm, steps, shield_mm, pinhole_dia_mm, fringe )
    # pr.trajectory_plot(ax, R_vec_mm)
    # print(R_vec_mm[-1])
    
    
    # mag = fields.Magenetic_field_Analitic(B_T[0], width_mm, depth_mm, k=0, fringe=fringe, sharp_edge=1, yoke=yoke)
    
    # fields.plot_magnetic_quiver(ax, mag, width_mm, depth_mm, height_mm, yoke,fringe, shield_mm)
    