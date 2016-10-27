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
x = np.loadtxt("x_field.txt")
y = np.loadtxt("y_field.txt")

X, Y = np.meshgrid(x,y)

# -- Load the data.
field = np.loadtxt("SalaminField.txt")


# -- Plot the field.
fig = plt.figure()
plt.pcolormesh(X,Y,field)
plt.gca().set_aspect('equal')
plt.colorbar()
plt.show()
