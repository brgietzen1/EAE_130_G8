#EAE130 Preliminary Sizing

# Standard Imports
import math
import numpy as np
import matplotlib.pyplot as plt


# 4.3 Drag Polar Estimates
# Functions
def find_wing_loading(rho, V, Cl):
    '''
    rho:    density    
    V:      velocity            [ft/s]
    Cl:     lift coefficient    [-]
    W_L :   wing loading        [lbs/ft^2]
    '''
    W_L = 0.5 * rho * V**2 * Cl
    return W_L

def drag_polar(W_L, C_f, AR, e, c, d, W_0, Cl, delta_Cd0=0.0):
    '''
    Calculate Cd for a given Cl, with configuration-specific delta_Cd0 added to Cd₀.
    
    Parameters:
        W_L:        Wing loading (lbf/ft²)
        C_f:        Skin friction coefficient
        AR:         Aspect ratio
        e:          Span efficiency factor
        c, d:       Raymer table constants
        W_0:        Gross weight (lbf)
        Cl:         Lift coefficient
        delta_Cd0:  Additional zero-lift drag coefficient for the configuration

    Returns:
        Cd:         Drag coefficient
    '''
    # Calculate wing area (S)
    S = W_0 / W_L
    
    # Calculate wetted area (S_wet)
    S_wet = 10 ** c * W_0 ** d
    
    # Calculate skin friction drag (f)
    f = C_f * S_wet
    
    # Base zero-lift drag coefficient (Cd₀)
    Cd0 = f / S
    
    # Total Cd₀ (base + configuration-specific delta_Cd0)
    Cd0_total = Cd0 + delta_Cd0
    
    # Total drag coefficient (Cd)
    Cd = Cd0_total + (Cl ** 2) / (math.pi * e * AR)
    
    return Cd

# Range of Cl values
def Cl_array(Cl_max):
    '''Generate an array of Cl values from -Cl_max to Cl_max'''
    return np.linspace(-Cl_max, Cl_max, 100)

# Create a dictionary for the different configuration parameters
configurations = [
    {
        "name": "Clean",
        "e": 0.85,
        "delta_Cd0": 0.0,  # No additional drag for clean configuration
        "Cl_max": 1.3, # Lift coefficient for clean
    },
    {
        "name": "Takeoff Flap Gear Up",
        "e": 0.8,
        "delta_Cd0": 0.01 ,
        "Cl_max": 1.5,
    },
    {
        "name": "Takeoff Flap Gear Down",
        "e": 0.8,
        "delta_Cd0": 0.025,
        "Cl_max": 1.5,
    },
    {
        "name": "Landing Flap Gear Up",
        "e": 0.75,
        "delta_Cd0": 0.055,
        "Cl_max": 1.7,
    },
    {
        "name": "Landing Flap Gear Down",
        "e": 0.75,
        "delta_Cd0": 0.055 + 0.015,
        "Cl_max": 1.7,
    }
]

# Parameters
# Fixed parameters (shared across all configurations)
W_L = 25       # Wing loading (lbf/ft²)
W_0 = 5000     # Gross weight (lbf)
C_f = 0.005    # Skin friction coefficient
AR = 8         # Aspect ratio (TBD)
c = 1.0447        # Roskam constant pg 122(134)
d = 0.5326        # Roskam constant pg 122(134)

# Plot setup
plt.figure(figsize=(10, 6))

#Cl_values = np.linspace(-1.3,1.3,100)
# Loop through each configuration
for config in configurations:
    name = config["name"]
    e = config["e"]
    delta_Cd0 = config["delta_Cd0"]
    Cl_max = config["Cl_max"]
    Cl_values = Cl_array(Cl_max)
    # Calculate Cd for all Cl values
    Cd_values = [drag_polar(W_L, C_f, AR, e, c, d, W_0, Cl, delta_Cd0) for Cl in Cl_values]
    #Cd_values = [drag_polar(W_L, C_f, AR, e, c, d, W_0, Cl, delta_Cd0)]
    # Plot the drag polar for this configuration
    plt.plot(Cd_values, Cl_values, label=name)

plt.xlabel("Drag Coefficient ($C_d$)")
plt.ylabel("Lift Coefficient ($C_l$)")
plt.title("Drag Polar for Different Aircraft Configurations")
plt.grid(True)
plt.legend()
plt.show()

# Constraint graph

# Stall speed
def stall_speed(rho, V_stall, Cl_max):
    '''find stall speed as a constant variables on the WS TW graph
    Parameters:
    rho: density
    V_stall: stall velocity
    Cl_max: max lift coefficient
    '''
    WS = 0.5 * rho * V_stall ** 2 * Cl_max

    return WS