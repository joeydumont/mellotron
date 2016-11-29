# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                    <joey.dumont@gmail.com>        #
# Created:      Oct. 19th, 2016                                               #
# Description:  We plot the Salamin fields computed in MELLOTRON.             #
# Dependencies: - NumPy                                                       #
#               - Matplotlib                                                  #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------ MAIN FUNCTION ------------------------------ #
# -- Define the mesh.
x_salamin = np.loadtxt("x_field_salamin.txt")/1e-6
y_salamin = np.loadtxt("y_field_salamin.txt")/1e-6

X_salamin, Y_salamin = np.meshgrid(x_salamin,y_salamin)

# -- Load the data.
field = np.loadtxt("SalaminField.txt")


# -- Plot the field.
fig = plt.figure()
plt.pcolormesh(X_salamin,Y_salamin,field)
plt.gca().set_aspect('equal')
plt.colorbar()


x_qgauss = np.loadtxt("x_field_qgauss.txt")/1e-6
y_qgauss = np.loadtxt("y_field_qgauss.txt")/1e-6
X_qgauss, Y_qgauss = np.meshgrid(x_qgauss, y_qgauss)
field = np.loadtxt("QuasiGaussianField.txt")

# -- Plot the field.
fig = plt.figure()
plt.pcolormesh(X_qgauss,Y_qgauss,field)
plt.gca().set_aspect('equal')
plt.colorbar()

t_data = np.loadtxt("t_data_qgauss.txt")
t_field = np.loadtxt("t_field_qgauss.txt")

fig = plt.figure()
plt.plot(t_data,t_field)

plt.show()