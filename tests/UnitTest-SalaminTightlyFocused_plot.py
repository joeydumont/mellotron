# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                    <joey.dumont@gmail.com>        #
# Created:      Oct. 19th, 2016                                               #
# Description:  We plot the Salamin fields computed in MELLOTRON.             #
# Dependencies: - NumPy                                                       #
#               - Matplotlib                                                  #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
import  matplotlib as mpl
mpl.use('pgf')
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

# ------------------------------ Configuration ------------------------------ #
#-- We reset the LaTeX parameters to enable XeLaTeX.
pgf_with_pdflatex = {
   "font.family": "serif", # use serif/main font for text elements
   "text.usetex": True,    # use inline math for ticks
   "pgf.rcfonts": False,   # don't setup fonts from rc parameters
   "pgf.preamble": [
        r"\usepackage{siunitx}",
        r"\usepackage{mathspec}",
        r"\usepackage[charter]{mathdesign}",
        r"\usepackage{fontspec}",
        r"\setmathfont{Fira Sans}",
        #r"\setmainfont{Oswald}",
        ]
}
mpl.rcParams.update(pgf_with_pdflatex)
mpl.rcParams['font.size'] = 10

# ------------------------------ MAIN FUNCTION ------------------------------ #
# -- Define the mesh.
x_salamin = np.loadtxt("x_field_salamin.txt")/800e-9
y_salamin = np.loadtxt("y_field_salamin.txt")/800e-9

X_salamin, Y_salamin = np.meshgrid(x_salamin,y_salamin)

# -- Load the data.
Ex_salamin = np.transpose(np.loadtxt("SalaminField_Ex.txt"))
Ey_salamin = np.transpose(np.loadtxt("SalaminField_Ey.txt"))
Ez_salamin = np.transpose(np.loadtxt("SalaminField_Ez.txt"))
Bx_salamin = np.transpose(np.loadtxt("SalaminField_Bx.txt"))
By_salamin = np.transpose(np.loadtxt("SalaminField_By.txt"))
Bz_salamin = np.transpose(np.loadtxt("SalaminField_Bz.txt"))

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.3)
plotOptions = {'rasterized':True, 'cmap': 'inferno'}

# ----------- Ex ------------- #
axEx = plt.subplot2grid((2,3), (0,0))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ex_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_x$")
axEx.set_ylabel(r"$y/\lambda$")


divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ey ------------- #
axEx = plt.subplot2grid((2,3), (0,1))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ey_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_y$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ez ------------- #
axEx = plt.subplot2grid((2,3), (0,2))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ez_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_z$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bx ------------- #
axEx = plt.subplot2grid((2,3), (1,0))
im   = plt.pcolormesh(X_salamin,Y_salamin,Bx_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_x$")
axEx.set_ylabel(r"$y/\lambda$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- By ------------- #
axEx = plt.subplot2grid((2,3), (1,1))
im   = plt.pcolormesh(X_salamin,Y_salamin,By_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_y$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bz ------------- #
axEx = plt.subplot2grid((2,3), (1,2))
im   = plt.pcolormesh(X_salamin,Y_salamin,Bz_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_z$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

plt.savefig("SalaminComponents.pdf", bbox_inches='tight', dpi=500)

# -------------------------------- E-DIPOLES -------------------------------- #
x_qgauss = np.loadtxt("x_field_qgauss.txt")/800.0e-9
y_qgauss = np.loadtxt("y_field_qgauss.txt")/800.0e-9
X_qgauss, Y_qgauss = np.meshgrid(x_qgauss, y_qgauss)

# -- Load the data.
Ex_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_Ex.txt"))
Ey_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_Ey.txt"))
Ez_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_Ez.txt"))
Bx_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_Bx.txt"))
By_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_By.txt"))
Bz_qgauss = np.transpose(np.loadtxt("QuasiGaussianField_Bz.txt"))

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.5)
plotOptions = {'rasterized':True, 'cmap': 'inferno'}

# ----------- Ex ------------- #
axEx = plt.subplot2grid((2,3), (0,0))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,Ex_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_x$")
axEx.set_ylabel(r"$y/\lambda$")


divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ey ------------- #
axEx = plt.subplot2grid((2,3), (0,1))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,Ey_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_y$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ez ------------- #
axEx = plt.subplot2grid((2,3), (0,2))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,Ez_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_z$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bx ------------- #
axEx = plt.subplot2grid((2,3), (1,0))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,Bx_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_x$")
axEx.set_ylabel(r"$y/\lambda$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- By ------------- #
axEx = plt.subplot2grid((2,3), (1,1))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,By_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_y$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bz ------------- #
axEx = plt.subplot2grid((2,3), (1,2))
im   = plt.pcolormesh(X_qgauss,Y_qgauss,Bz_qgauss, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_z$")
axEx.set_xlabel(r"$x/\lambda$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

plt.savefig("QuasiGaussianComponents.pdf", bbox_inches='tight', dpi=500)