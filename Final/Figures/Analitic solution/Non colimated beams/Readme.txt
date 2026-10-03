 me_kg = 9.109*1e-31
    R0_mm  = np.array([0,-15,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = -5e7
    yoke = 0
    shield_mm = 5
    pinhole_dia_mm = 3
    fringe = 1
    CFL = 0.05
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =    "Analitic field" #  "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**8 + 1;
    Ny_p = 2**8 + 1;
    

    # spec = Spectrometer_Budy(shield_mm = shield_mm, yoke=yoke)
    # fig,ax = spec._draw_spec()
    
    experiment = Experiment(
                            height_mm = height_mm, width_mm = width_mm, depth_mm = depth_mm, shield_mm = shield_mm, pinhole_dia_mm = pinhole_dia_mm, yoke = yoke, # Spectrometere and simulation border
                            R0_mm= R0_mm, q_eng_MeV = e_eng_MeV, m_kg = me_kg, q_C = e_C, # Particle parameters
                            Bx0_T = Bx0_T, Ex0_Vm = Ex0_Vm, fringe = fringe, sharp_edge = sharp_edge, # Fildes parameters
                            CFL = CFL, N_steps = N_steps, solution = solution, # Simulation resolutions
                            Nx_p = Nx_p, Ny_p = Ny_p # Potential solving
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
        integrator.psi._show_potential()
        integrator.mag._show_field()
        integrator.phi._show_potential()
        integrator.Ele._show_field()
        
    




Figure 1:

        integrators_name = ["RK4_Herm"]#["Boris","RK2","RK4","Euler"]
        params["fringe"] = 1
        params["Ex0_Vm"] = 0
        params["R0_mm"] = np.array([0,-3,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        
        e_engs_MeV_scalar = 8
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Exponential"#"Uniform" #"Gaussian" #"Exponential"


Figure 2:

integrators_name = ["RK4_Herm"]#["Boris","RK2","RK4","Euler"]
        params["fringe"] = 1
        params["Ex0_Vm"] = -1e7
        params["R0_mm"] = np.array([0,-3,0])
        params = experiment._update_params(params)
        integrator._update_params(params)
        
        e_engs_MeV_scalar = 8
        N_steps_vec = 5e2#,4e3,5e3,1e4,3e4,5e4]
        number_of_experiments = 1000
        angular_distribution = "Exponential"#"Uniform" #"Gaussian" #"Exponential"