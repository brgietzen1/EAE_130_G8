import pyxfoil
import numpy as np
import matplotlib.pyplot as plt

# Getting Cp data for both inviscid and viscous cases
naca = True
alf = np.linspace(0,21)
foil = '0016'
Re = 5e6
pyxfoil.GetPolar(foil, naca, alf, Re)

# Extracting x and y coordinates of NACA 0014
filename = 'Data/naca0016/naca0016.dat'
x, y = 4.63* np.loadtxt(filename, dtype=float, unpack=True, skiprows = 1)

# plot the geometry and the panels
width = 10
plt.figure(figsize=(width, width))
plt.xlabel('x (ft)', fontsize=16)
plt.ylabel('y (ft)', fontsize=16)
plt.plot(x, y, color='b', linestyle='-', linewidth=2)
plt.title('NACA 0016 Geometry', fontsize = 16)
plt.axis('equal')
plt.grid(True)
plt.show()
