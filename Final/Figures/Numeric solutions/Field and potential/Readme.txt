 
    me_kg = 9.109*1e-31
    R0_mm  = np.array([0,0,0])
    e_C = -1.602*1e-19
    e_eng_MeV = np.array([0,10,0])
    height_mm =26; width_mm = 12.5; depth_mm   = 50.8
    Bx0_T = 0.5
    Ex0_Vm = 1e8
    yoke = 1
    shield_mm = 0
    pinhole_dia_mm = 3
    fringe = 0
    CFL = 0.05
    N_steps = 5*1e3
    pinhole_dia_mm = 3
    solution =    "Analitic field" #  "Numeric Field" # "Analitic field" #  
    rand = "uni" # exp #gauss # "None"
    sharp_edge = 1
    
    Nx_p = 2**8 + 1;
    Ny_p = 2**8 + 1;

    if False:#False:
        params["solution"] = "Analitic field"
    
        integrator = Integrators(params)
    
    if True:#True:
        params["solution"] = "Numeric Field"
        integrator = Integrators(params)
        integrator.psi._show_potential()
        integrator.mag._show_field()
        integrator.phi._show_potential()
        integrator.Ele._show_field()
        


residual value: 4.028615236961741e-07
 number of iterations: 130
residual value: 4.028615236961741e-07
 number of iterations: 130