# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                    <joey.dumont@gmail.com>        #
# Created:      Oct. 19th, 2016                                               #
# Description:  We plot the Salamin fields computed in MELLOTRON.             #
# Dependencies: - NumPy                                                       #
#               - Matplotlib                                                  #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
from mpl_toolkits.mplot3d import Axes3D
import  matplotlib as mpl
mpl.use('pgf')
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PolyCollection

import numpy as np
import scipy.signal as sig
import scipy.constants as cst
import numpy.fft as fft
import math
import matplotlib.animation as animation
import shlex
import subprocess

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
x_salamin = np.loadtxt("x_field_salamin.txt")*1e6
y_salamin = np.loadtxt("y_field_salamin.txt")*1e6

X_salamin, Y_salamin = np.meshgrid(x_salamin,y_salamin)

# -- Load the data.
Ex_salamin = np.transpose(np.loadtxt("SalaminField_Ex.txt"))**2
Ey_salamin = np.transpose(np.loadtxt("SalaminField_Ey.txt"))**2
Ez_salamin = np.transpose(np.loadtxt("SalaminField_Ez.txt"))**2
Bx_salamin = np.transpose(np.loadtxt("SalaminField_Bx.txt"))**2
By_salamin = np.transpose(np.loadtxt("SalaminField_By.txt"))**2
Bz_salamin = np.transpose(np.loadtxt("SalaminField_Bz.txt"))**2

max_field = np.amax([Ex_salamin, Ey_salamin, Ez_salamin, Bx_salamin, By_salamin, Bz_salamin])

Ex_salamin /= max_field
Ey_salamin /= max_field
Ez_salamin /= max_field
Bx_salamin /= max_field
By_salamin /= max_field
Bz_salamin /= max_field

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.3)
plotOptions = {'rasterized':True, 'cmap': 'viridis'}

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

# ------------------------------ Waterfall Plot ----------------------------- #
salamin_time_times = np.loadtxt("SalaminTimeIe.txt")
size_time       = salamin_time_times.size
Ex_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))
Ey_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))
Ez_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))
Bx_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))
By_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))
Bz_salamin_time = np.empty((Ex_salamin.shape[0], Ex_salamin.shape[1], size_time))

# -- Load the data.
for i in range(size_time):
  Ex_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_Ex_t{}.txt".format(i)))
  Ey_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_Ey_t{}.txt".format(i)))
  Ez_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_Ez_t{}.txt".format(i)))
  Bx_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_Bx_t{}.txt".format(i)))
  By_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_By_t{}.txt".format(i)))
  Bz_salamin_time[:,:,i] = np.transpose(np.loadtxt("SalaminField_Bz_t{}.txt".format(i)))

# -- Prepare the waterfall plot.
# -- Prepare a LogNorm for colours and alpha.
x_cut_time      = x_salamin
field_xcut_time = Ex_salamin_time[Ex_salamin_time.shape[0]//2  ,:,:]**2 \
                  + Ey_salamin_time[Ey_salamin_time.shape[0]//2,:,:]**2 \
                  + Ez_salamin_time[Ez_salamin_time.shape[0]//2,:,:]**2
longfield_xcut  = Ez_salamin_time[Ez_salamin_time.shape[0]//2,:,:]**2
field_xcut_time /= np.amax(field_xcut_time)
longfield_xcut  /= np.amax(longfield_xcut)
time            = salamin_time_times/1e-15
norm = mpl.colors.PowerNorm(gamma=0.1,vmin=0.0, vmax=np.amax(np.abs(time)))

fig = plt.figure(figsize=(4,4))
ax  = fig.add_subplot(111, projection='3d')
ax.view_init(20,110)
verts = []
vertsl = []
max_field = np.empty((size_time))
max_fieldl= np.empty((size_time))

for i in range(size_time):
    xs = np.concatenate([[x_cut_time[0]], x_cut_time[:,], [x_cut_time[-1]]])
    ys = np.concatenate([[0],field_xcut_time[:,i],[0]])
    ysl= np.concatenate([[0],longfield_xcut[:,i],[0]])
    verts.append(list(zip(xs,ys)))
    vertsl.append(list(zip(xs,ysl)))
    max_field[i] = np.amax(field_xcut_time[:,i])
    max_fieldl[i]= np.amax(longfield_xcut[:,i])
    ax.plot(x_cut_time[:],field_xcut_time[:,i], zs=time[-1], zdir='y', zorder=-1, color='C0', alpha=1.0-norm(np.abs(time[i])))

poly = PolyCollection(verts, rasterized=True, facecolor=None, edgecolor='k', lw=0.7)
polyl= PolyCollection(vertsl, rasterized=True, facecolor=None, edgecolor='k', lw=0.7)
poly.set_alpha(0.25)
polyl.set_alpha(0.25)
#ax.scatter(time, max_field, zs=x_cut_time[0,0], zdir='x', zorder=-1)
#ax.plot(time, EnvelopeFunction(time,*popt), zs=x_cut_time[0,0], zdir='x', zorder=-1)
#ax.plot(x_cut_time[:,time_idx], field_xcut_time[:,time_idx], zs=time[-1], zdir='y', zorder=-1)
ax.add_collection3d(poly, zs=time, zdir='y')
#ax.add_collection3d(polyl,zs=time,zdir='y')
ax.ticklabel_format(style='sci', axis='z',scilimits=(0,0))
ax.set_xlim3d(x_cut_time.min(), x_cut_time.max())
ax.set_xlabel(r'$x$ [$\mu m$]')
ax.set_ylim3d(np.amin(time), np.amax(time))
ax.set_ylabel('Time (fs)')
ax.invert_yaxis()
ax.invert_xaxis()
ax.set_zlim3d(0.0, field_xcut_time.max())
ax.set_zlabel('Amplitude', rotation=90)
plt.savefig("SalaminElectricIntensityTimeWaterfall.pdf",dpi=500)
plt.close()

plt.figure()
plt.plot(time, max_field)
plt.plot(time, max_fieldl)
plt.savefig("SalaminIevsEz.pdf", dpi=500)
plt.close()

# -- Further processing of the waterfall figure.
cmd = 'pdfcrop --margins "35 0 0 0" SalaminElectricIntensityTimeWaterfall.pdf SalaminElectricIntensityTimeWaterfall-cropped.pdf'
proc = subprocess.call(shlex.split(cmd))
print("Waterfall plot done.")
# ------------------------------ FFT of Salamin ----------------------------- #
Ex_salamin_freq = fft.rfft(Ex_salamin_time)
Ey_salamin_freq = fft.rfft(Ey_salamin_time)
Ez_salamin_freq = fft.rfft(Ez_salamin_time)
Bx_salamin_freq = fft.rfft(Bx_salamin_time)
By_salamin_freq = fft.rfft(By_salamin_time)
Bz_salamin_freq = fft.rfft(Bz_salamin_time)

freqs       = fft.rfftfreq(size_time, d=salamin_time_times[1]-salamin_time_times[0])
wavelengths = np.pi*cst.c/freqs

print(wavelengths*1e9)

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1, array[idx-1]
    else:
        return idx, array[idx]
idx, actual_freq = find_nearest(wavelengths[::-1]*1e9, 800.0)
print("Actual freq is {}".format(actual_freq))

plt.figure()
plt.pcolormesh(X_salamin, Y_salamin, np.abs(Ex_salamin_time[:,:,size_time//2])**2, **plotOptions)
plt.gca().set_aspect('equal')
plt.savefig("TestTime.pdf", bbox_inches='tight', dpi=500)
plt.close()

# -- Plot intensity of Ex.
intensity        = np.abs(Ex_salamin_freq)**2
maxIntensity_idx = np.argmax(intensity)
maxIntensity_ind = np.unravel_index(maxIntensity_idx, intensity.shape)
maxIntensity     = intensity[maxIntensity_ind]

plt.figure()
im = plt.pcolormesh(X_salamin, Y_salamin, intensity[:,:,maxIntensity_ind[2]]/maxIntensity, vmin=0.0, vmax=1.0,**plotOptions)
plt.gca().set_aspect('equal')
plt.gca().text(0.95,0.95, "freq = {:5e}".format(2*np.pi*freqs[maxIntensity_ind[2]]))
plt.savefig("TestFFT.pdf", bbox_inches='tight', dpi=500)

# -------------------------- Frequency Components -------------------------- #
Ex_intensity = np.abs(Ex_salamin_freq)**2
Ey_intensity = np.abs(Ey_salamin_freq)**2
Ez_intensity = np.abs(Ez_salamin_freq)**2
Bx_intensity = np.abs(Bx_salamin_freq)**2
By_intensity = np.abs(By_salamin_freq)**2
Bz_intensity = np.abs(Bz_salamin_freq)**2

Ex_intensity_ind = np.unravel_index(np.argmax(Ex_intensity), Ex_intensity.shape)
Ey_intensity_ind = np.unravel_index(np.argmax(Ey_intensity), Ex_intensity.shape)
Ez_intensity_ind = np.unravel_index(np.argmax(Ex_intensity), Ex_intensity.shape)
Bx_intensity_ind = np.unravel_index(np.argmax(Bx_intensity), Ex_intensity.shape)
By_intensity_ind = np.unravel_index(np.argmax(By_intensity), Ex_intensity.shape)
Bz_intensity_ind = np.unravel_index(np.argmax(Bz_intensity), Ex_intensity.shape)

listofInd = [Ex_intensity_ind, Ey_intensity_ind, Ez_intensity_ind, Bx_intensity_ind, By_intensity_ind, Bz_intensity_ind]
listOfMax = [Ex_intensity[Ex_intensity_ind], Ey_intensity[Ey_intensity_ind], Ez_intensity[Ez_intensity_ind], \
                    Bx_intensity[Bx_intensity_ind], By_intensity[By_intensity_ind], Bz_intensity[Bz_intensity_ind]]
maxInd  = np.argmax(listOfMax)

Ex_intensity /= listOfMax[maxInd]
Ey_intensity /= listOfMax[maxInd]
Ez_intensity /= listOfMax[maxInd]
Bx_intensity /= listOfMax[maxInd]
By_intensity /= listOfMax[maxInd]
Bz_intensity /= listOfMax[maxInd]

Ex_intensity_plot = Ex_intensity[:,:,listofInd[maxInd][2]]
Ey_intensity_plot = Ey_intensity[:,:,listofInd[maxInd][2]]
Ez_intensity_plot = Ez_intensity[:,:,listofInd[maxInd][2]]
Bx_intensity_plot = Bx_intensity[:,:,listofInd[maxInd][2]]
By_intensity_plot = By_intensity[:,:,listofInd[maxInd][2]]
Bz_intensity_plot = Bz_intensity[:,:,listofInd[maxInd][2]]

# -- Plot the field.
fig  = plt.figure(figsize=(7,4))
fig.subplots_adjust(wspace=0.5,hspace=0.3)
plotOptions = {'rasterized':True, 'cmap': 'viridis'}

# ----------- Ex ------------- #
axEx = plt.subplot2grid((2,3), (0,0))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ex_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_x$")
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")


divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ey ------------- #
axEx = plt.subplot2grid((2,3), (0,1))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ey_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_y$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Ez ------------- #
axEx = plt.subplot2grid((2,3), (0,2))
im   = plt.pcolormesh(X_salamin,Y_salamin,Ez_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$E_z$")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bx ------------- #
axEx = plt.subplot2grid((2,3), (1,0))
im   = plt.pcolormesh(X_salamin,Y_salamin,Bx_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_x$")
axEx.set_ylabel(r"$y\,[\mu\si{\micro\metre}]")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- By ------------- #
axEx = plt.subplot2grid((2,3), (1,1))
im   = plt.pcolormesh(X_salamin,Y_salamin,By_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_y$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

# ----------- Bz ------------- #
axEx = plt.subplot2grid((2,3), (1,2))
im   = plt.pcolormesh(X_salamin,Y_salamin,Bz_intensity_plot, **plotOptions)
axEx.set_aspect('equal')
axEx.set_title(r"$B_z$")
axEx.set_xlabel(r"$x\,[\mu\si{\metre}]")

divider = make_axes_locatable(axEx)
cax     = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax)

plt.savefig("SalaminComponentsFreq.pdf", bbox_inches='tight', dpi=500)

# -- Animation of E_x(k).
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

freqFigAnimation = plt.figure()
ims = []

for i in range(freqs.size):
  im = plt.pcolormesh(X_salamin, Y_salamin, np.abs(Ex_salamin_freq[:,:,i])**2/maxIntensity, vmin=0.0, vmax=1.0,**plotOptions)
  t  = plt.gca().text(0.95,0.95, "wavelength = {:.5e}".format(wavelengths[i]*1e9), transform=plt.gca().transAxes)
  f  = plt.gca().text(0.95,0.85, "freq       = {:.5e}".format(freqs[i]),           transform=plt.gca().transAxes)
  ims.append([im,t,f])
  #plt.savefig("TestFFT.pdf", bbox_inches='tight',dpi=500)

plt.gca().set_aspect('equal')
im_ani = animation.ArtistAnimation(freqFigAnimation, ims, interval=500, blit=True, repeat=False)
im_ani.save("im.mp4", writer=writer)

plt.close()

# ----------------------------- Salamin in Time ----------------------------- #
salamin_time_intensity     = np.loadtxt("SalaminTimeIe_time.txt")
salamin_time_times = np.loadtxt("SalaminTimeIe.txt")/1e-15

# -- Compute the Hilbert transform of this.
salamin_envelope = np.abs(sig.hilbert(salamin_time_intensity))
fig = plt.figure()
ax  = fig.add_subplot(111)
im  = plt.plot(salamin_time_times, salamin_time_intensity, rasterized=True)
plt.plot(salamin_time_times, salamin_envelope, 'k--')

ax.set_xlabel(r"Time [fs]")
ax.set_ylabel(r"Intensity [$\si{\watt\per\cm\squared}$]")

plt.savefig("SalaminTimeIntensity.pdf", bbox_inches='tight')

# # ------------------------- Quasi-Gaussian e-dipoles ------------------------ #
# x_qgauss = np.loadtxt("x_field_qgauss.txt")/1e-6
# y_qgauss = np.loadtxt("y_field_qgauss.txt")/1e-6
# X_qgauss, Y_qgauss = np.meshgrid(x_qgauss, y_qgauss)
# field = np.loadtxt("QuasiGaussianField.txt")
# fieldQED=np.loadtxt("QuasiGaussianField_qed.txt")

# # -- Plot the field.
# fig = plt.figure()
# plt.pcolormesh(X_qgauss,Y_qgauss,field)
# plt.gca().set_aspect('equal')
# plt.colorbar()

# fig = plt.figure()
# plt.pcolormesh(X_qgauss,Y_qgauss,(field-fieldQED)/field, **plotOptions)
# plt.gca().set_aspect('equal')
# plt.colorbar()

# t_data = np.loadtxt("t_data_qgauss.txt")
# t_field = np.loadtxt("t_field_qgauss.txt")
# t_field_qed= np.loadtxt("t_field_qgauss_qed.txt")

# fig = plt.figure()
# plt.plot(t_data,t_field)
# plt.plot(t_data,t_field_qed)


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
