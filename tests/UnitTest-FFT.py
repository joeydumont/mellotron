# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                    <joey.dumont@gmail.com>        #
# Created:      Jul. 24th, 2017                                               #
# Description:  We test different FFT implementations.                        #
# Dependencies: - NumPy                                                       #
#               - Matplotlib                                                  #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
import  matplotlib as mpl
#mpl.use('pgf')
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import scipy.fftpack as fft_sp

# ------------------------------ Configuration ------------------------------ #
#-- We reset the LaTeX parameters to enable XeLaTeX.
# pgf_with_pdflatex = {
#    "font.family": "serif", # use serif/main font for text elements
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.preamble": [
#         r"\usepackage{siunitx}",
#         r"\usepackage{mathspec}",
#         r"\usepackage[charter]{mathdesign}",
#         r"\usepackage{fontspec}",
#         r"\setmathfont{Fira Sans}",
#         #r"\setmainfont{Oswald}",
#         ]
# }
# mpl.rcParams.update(pgf_with_pdflatex)
# mpl.rcParams['font.size'] = 10

# ------------------------------ MAIN FUNCTION ------------------------------ #

# Samples
size = 500
x    = np.linspace(-10,10,size)
dt   = x[1]-x[0]
y    = np.exp(-x**2)

# Compute FFT of samples with NumPy.
x_fft_numpy     = 2*np.pi*np.fft.rfftfreq(size,d=dt)
y_fft_numpy     = dt*np.fft.rfft(y)/np.sqrt(2*np.pi)
y_fft_numpy_ana = np.exp(-(x_fft_numpy)**2/4)/np.sqrt(2)

# Compute FFT of samples with SciPy.
x_fft_scipy = 2*np.pi*fft_sp.rfftfreq(size, d=dt)
y_fft_scipy = dt*fft_sp.rfft(y)/np.sqrt(2*np.pi)

x_fft_scipy_spliced = np.concatenate([[x_fft_scipy[0]], x_fft_scipy[1:-1:2]])
y_fft_scipy_spliced = np.concatenate([[y_fft_scipy[0]], y_fft_scipy[1:-1:2]])+1j*np.concatenate([[0], y_fft_scipy[2:-1:2]])
y_fft_scipy_ana     = np.exp(-x_fft_scipy_spliced**2/4)/np.sqrt(2)

# Compare plots.
plt.figure()
plt.plot(x_fft_numpy,         np.abs(y_fft_numpy), label="NumPy")
plt.plot(x_fft_numpy,         y_fft_numpy_ana,     label="NumPy - analytical")
plt.legend(loc=0)

plt.figure()
plt.plot(x_fft_scipy_spliced, np.abs(y_fft_scipy_spliced), label="SciPy")
plt.plot(x_fft_scipy_spliced, y_fft_scipy_ana,     label="SciPy - analytical")
plt.legend(loc=0)
plt.show()

