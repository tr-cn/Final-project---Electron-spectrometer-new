import numpy as np

import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'mathtext.fontset': 'stix'  # Matches Times New Roman math styles seamlessly
})

from matplotlib.ticker import MaxNLocator, FormatStrFormatter


def trajectory_plot (ax, R_vec_mm):
    R_vec_mm = np.array(R_vec_mm)
    x_vec_mm = R_vec_mm[:,0]; z_vec_mm = R_vec_mm[:,1]; y_vec_mm = R_vec_mm[:,2];
    
    ax.plot(x_vec_mm, z_vec_mm, y_vec_mm)
    
    return 

