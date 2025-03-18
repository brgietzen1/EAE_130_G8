import pyxfoil
import numpy as np
import matplotlib.pyplot as plt

Re = 5e6
alpha_range = np.arange(0, 22, 0.25)  # Angles of attack from 0 to 21 degrees

# Run XFOIL to generate the polar
obj = pyxfoil.GetPolar(foil='23018', naca=True, alfs=alpha_range, Re=Re, SaveCP=False)

# Load the generated polar file (assuming it is stored in 'Polars/' directory)
polar_filename = f'Data/naca23018/naca23018_polar_Re5.00e+06a0.0-21.8.dat'
polar_data = np.loadtxt(polar_filename, skiprows=12)

# Extract angle of attack (alpha), lift coefficient (Cl), and drag coefficient (Cd)
alpha = polar_data[:, 0]
Cl = polar_data[:, 1]
Cd = polar_data[:, 2]

# Find and print the maximum Cl value
max_Cl = np.max(Cl)
max_Cl_alpha = alpha[np.argmax(Cl)]
print(f"Maximum Cl: {max_Cl:.4f} at α = {max_Cl_alpha:.2f}°")

# Plot Cl vs. Alpha (Lift Curve)
plt.figure(figsize=(8, 6))
plt.plot(alpha, Cl, marker='o', linestyle='-')
plt.xlabel("Angle of Attack (degrees)")
plt.ylabel("Lift Coefficient $C_l$")
plt.title("Lift Curve for NACA 23018 (Re=5×10⁶)")
plt.grid(True)
plt.show()

