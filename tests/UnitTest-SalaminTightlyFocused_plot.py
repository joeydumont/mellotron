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
import scipy.signal as sig
import scipy.constants as cst

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
Ex_salamin = np.transpose(np.loadtxt("SalaminField_Ex.txt"))**2
Ey_salamin = np.transpose(np.loadtxt("SalaminField_Ey.txt"))**2
Ez_salamin = np.transpose(np.loadtxt("SalaminField_Ez.txt"))**2
Bx_salamin = np.transpose(np.loadtxt("SalaminField_Bx.txt"))**2
By_salamin = np.transpose(np.loadtxt("SalaminField_By.txt"))**2
Bz_salamin = np.transpose(np.loadtxt("SalaminField_Bz.txt"))**2

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.3)
plotOptions = {'rasterized':True, 'cmap': 'inferno'}

# ----------- Ex ------------- #
axEx = plt.subplot2grid((2,3), (0,0))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ex_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_x$")
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")


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
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- By ------------- #
axEx = plt.subplot2grid((2,3), (1,1))
im   = plt.pcolormesh(X_salamin,Y_salamin,By_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_y$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bz ------------- #
axEx = plt.subplot2grid((2,3), (1,2))
im   = plt.pcolormesh(X_salamin,Y_salamin,Bz_salamin, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_z$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

plt.savefig("SalaminComponents.pdf", bbox_inches='tight', dpi=500)

# ----------------------------- Salamin in Time ----------------------------- #
salamin_time_times     = np.loadtxt("SalaminTimeIe_time.txt")/1e-15
salamin_time_intensity = np.loadtxt("SalaminTimeIe.txt")

fig = plt.figure()
ax  = fig.add_subplot(111)
im  = plt.semilogy(salamin_time_times,salamin_time_intensity, '--o')
ax.set_xlabel(r"Time [fs]")
ax.set_ylabel(r"Intensity [$\si{\watt\per\cm\squared}$]")

plt.savefig("SalaminTimeIntensity.pdf", bbox_inches='tight')

x_qgauss = np.loadtxt("x_field_qgauss.txt")/1e-6
y_qgauss = np.loadtxt("y_field_qgauss.txt")/1e-6
X_qgauss, Y_qgauss = np.meshgrid(x_qgauss, y_qgauss)
field = np.loadtxt("QuasiGaussianField.txt")
fieldQED=np.loadtxt("QuasiGaussianField_qed.txt")

# ------------------------- Quasi-Gaussian e-dipoles ------------------------ #
# -- Plot the field.
fig = plt.figure()
plt.pcolormesh(X_qgauss,Y_qgauss,field)
plt.gca().set_aspect('equal')
plt.colorbar()

fig = plt.figure()
plt.pcolormesh(X_qgauss,Y_qgauss,(field-fieldQED)/field)
plt.gca().set_aspect('equal')
plt.colorbar()

t_data = np.loadtxt("t_data_qgauss.txt")
t_field = np.loadtxt("t_field_qgauss.txt")
t_field_qed= np.loadtxt("t_field_qgauss_qed.txt")

fig = plt.figure()
plt.plot(t_data,t_field)
plt.plot(t_data,t_field_qed)


# ------------------------- StrattoCalculator Linear ------------------------ #
# -- Define the mesh.
lambda_c = 800.0e-9
omega_0  = 2*np.pi*cst.c/lambda_c
x_stratto = np.loadtxt("x_field_stratto.txt")/lambda_c
y_stratto = np.loadtxt("y_field_stratto.txt")/lambda_c

X_stratto, Y_stratto = np.meshgrid(x_stratto,y_stratto)

# -- Load the data.
intensity_prefac = 0.5*cst.c*cst.epsilon_0*(cst.m_e*omega_0*cst.c/cst.e)**2
Ex_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_Ex.txt"))**2
Ey_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_Ey.txt"))**2
Ez_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_Ez.txt"))**2
Bx_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_Bx.txt"))**2
By_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_By.txt"))**2
Bz_stratto = 1.0e-4*intensity_prefac*np.transpose(np.loadtxt("StrattoField_Bz.txt"))**2

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.3)
plotOptions = {'rasterized':True, 'cmap': 'inferno'}

# ----------- Ex ------------- #
axEx = plt.subplot2grid((2,3), (0,0))
im   = plt.pcolormesh(X_stratto,Y_stratto,Ex_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_x$")
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")


divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ey ------------- #
axEx = plt.subplot2grid((2,3), (0,1))
im   = plt.pcolormesh(X_stratto,Y_stratto,Ey_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_y$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ez ------------- #
axEx = plt.subplot2grid((2,3), (0,2))
im   = plt.pcolormesh(X_stratto,Y_stratto,Ez_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_z$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bx ------------- #
axEx = plt.subplot2grid((2,3), (1,0))
im   = plt.pcolormesh(X_stratto,Y_stratto,Bx_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_x$")
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- By ------------- #
axEx = plt.subplot2grid((2,3), (1,1))
im   = plt.pcolormesh(X_stratto,Y_stratto,By_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_y$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bz ------------- #
axEx = plt.subplot2grid((2,3), (1,2))
im   = plt.pcolormesh(X_stratto,Y_stratto,Bz_stratto, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_z$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

plt.savefig("StrattoComponents.pdf", bbox_inches='tight', dpi=500)

# ----------------------------- StrattoLinear in Time ----------------------------- #
stratto_Ex_time     = np.loadtxt("StrattoField_Ex_time.txt")

fig = plt.figure()
ax  = fig.add_subplot(111)
im  = plt.plot(stratto_Ex_time[:,0],stratto_Ex_time[:,1], '--o')
#ax.set_xlabel(r"Time [fs]")
#ax.set_ylabel(r"Intensity [$\si{\watt\per\cm\squared}$]")

plt.savefig("Stratto_Ex_time.pdf", bbox_inches='tight')

x_qgauss = np.loadtxt("x_field_qgauss.txt")/1e-6
y_qgauss = np.loadtxt("y_field_qgauss.txt")/1e-6
X_qgauss, Y_qgauss = np.meshgrid(x_qgauss, y_qgauss)
field = np.loadtxt("QuasiGaussianField.txt")
fieldQED=np.loadtxt("QuasiGaussianField_qed.txt")

plt.show()