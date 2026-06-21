import matplotlib.pyplot as plt
import numpy as np

from IPython import get_ipython
plt.close('all'); get_ipython().run_line_magic('clear', ''); get_ipython().run_line_magic('reset', '-f');
from spectrometer.geometry import geometry


if __name__ == "__main__":
    geometry()
    