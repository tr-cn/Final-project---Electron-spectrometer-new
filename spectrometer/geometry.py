import numpy as np

import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'mathtext.fontset': 'stix'  # Matches Times New Roman math styles seamlessly
})

from matplotlib.ticker import MaxNLocator, FormatStrFormatter


def  geometry (height_mm =26, width_mm = 12.5, depth_mm   = 50.8, pinhole_dia_mm =3, shield_mm = 0):
    pinhole_radius_mm = pinhole_dia_mm/2
    fig = plt.figure()
    
    manager = plt.get_current_fig_manager()
    try:
        # Works perfectly on Spyder's default backend (Qt) and Windows/Linux/Mac
        manager.window.showMaximized()
    except AttributeError:
        try:
            # Fallback syntax for alternative backends (Tkinter)
            manager.resize(*manager.window.maxsize())
        except Exception:
            pass
    ax = fig.add_subplot(111, projection='3d')
    

   
    
    padding = abs(height_mm - width_mm)/2
    if height_mm-width_mm>=0:
        norm = np.sqrt(height_mm**2 + (width_mm + padding)**2 + depth_mm**2 )
        ax.set_box_aspect([(width_mm + padding)/norm, depth_mm/norm, height_mm/norm])
        ax.set_xlim([-width_mm/2-padding/2, width_mm/2 + padding/2])
    else:
        norm = np.sqrt((height_mm + padding)**2 + width_mm**2 + depth_mm**2 )
        ax.set_box_aspect([width_mm/norm, depth_mm/norm, (height_mm + padding)/norm])
        ax.set_zlim([-height_mm/2-padding/2, height_mm/2 + padding/2])
        
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)    

    # E-spec entrance 
    ax.plot([-width_mm/2, width_mm/2],   [0, 0], [-height_mm/2, -height_mm/2], color='black', linewidth=2)
    ax.plot([-width_mm/2, width_mm/2],   [0, 0], [height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([-width_mm/2, -width_mm/2], [0, 0], [-height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([width_mm/2, width_mm/2],   [0, 0], [-height_mm/2, height_mm/2],   color='black', linewidth=2)

    # E-spec end 
    ax.plot([-width_mm/2, width_mm/2],   [depth_mm, depth_mm], [-height_mm/2, -height_mm/2], color='black', linewidth=2)
    ax.plot([-width_mm/2, width_mm/2],   [depth_mm, depth_mm], [height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([-width_mm/2, -width_mm/2], [depth_mm, depth_mm], [-height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([width_mm/2, width_mm/2],   [depth_mm, depth_mm], [-height_mm/2, height_mm/2],   color='black', linewidth=2)

    # E-spec body 
    ax.plot([-width_mm/2, -width_mm/2], [0, depth_mm], [-height_mm/2, -height_mm/2], color='black', linewidth=2)
    ax.plot([-width_mm/2, -width_mm/2], [0, depth_mm], [height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([width_mm/2, width_mm/2],   [0, depth_mm], [height_mm/2, height_mm/2],   color='black', linewidth=2)
    ax.plot([width_mm/2, width_mm/2],   [0, depth_mm], [-height_mm/2, -height_mm/2], color='black', linewidth=2)

    # E-spec entry hole
    theta = np.linspace(0, 2*np.pi, 1001)
    ax.plot(pinhole_radius_mm * np.cos(theta), np.zeros_like(theta), pinhole_radius_mm * np.sin(theta), color='black', linewidth=2)
    
    if shield_mm>0:
        # Shiled face
        ax.plot([-width_mm/2, width_mm/2],   [-shield_mm, -shield_mm], [-height_mm/2, -height_mm/2], color='gray', linewidth=2)
        ax.plot([-width_mm/2, width_mm/2],   [-shield_mm, -shield_mm], [height_mm/2, height_mm/2],   color='gray', linewidth=2)
        ax.plot([-width_mm/2, -width_mm/2], [-shield_mm, -shield_mm], [-height_mm/2, height_mm/2],   color='gray', linewidth=2)
        ax.plot([width_mm/2, width_mm/2],   [-shield_mm, -shield_mm], [-height_mm/2, height_mm/2],   color='gray', linewidth=2)

  
        # Shiled body
        ax.plot([-width_mm/2, -width_mm/2],  [-shield_mm, 0], [-height_mm/2, -height_mm/2], color='gray', linewidth=2)
        ax.plot([-width_mm/2, -width_mm/2], [-shield_mm, 0], [height_mm/2, height_mm/2],   color='gray', linewidth=2)
        ax.plot([width_mm/2, width_mm/2],   [-shield_mm, 0], [height_mm/2, height_mm/2],   color='gray', linewidth=2)
        ax.plot([width_mm/2, width_mm/2],   [-shield_mm, 0], [-height_mm/2, -height_mm/2], color='gray', linewidth=2)

        #Shiled entry hole
        theta = np.linspace(0, 2*np.pi, 1001)
        ax.plot(pinhole_radius_mm * np.cos(theta), np.zeros_like(theta)-shield_mm, pinhole_radius_mm * np.sin(theta), color='gray', linewidth=2)




    lable_fontsize = 30
    tic_fontsize = 25
    padding_distance = 25
    ax.set_xlabel('Width (mm)',fontsize = lable_fontsize, labelpad=padding_distance); 
    ax.set_ylabel('Depth (mm)',fontsize = lable_fontsize,labelpad=padding_distance)
    ax.set_zlabel('Height (mm)',fontsize = lable_fontsize, labelpad=padding_distance)
    
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.zaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    ax.tick_params(axis='both', labelsize = tic_fontsize,  )

    ax.view_init(elev=15, azim=335)
    ax.dist = 9
    #plt.show()
    
    return fig, ax
