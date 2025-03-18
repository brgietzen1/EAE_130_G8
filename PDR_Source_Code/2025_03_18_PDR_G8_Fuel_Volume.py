import numpy as np

def CalcFuelVolume(S, b, gamma, t_ct, t_cr):
    """
    Calculates the estimated fuel volume capacity within a wing.
    
    Parameters:
    S (float): Wing area in square feet (ft²)
    b (float): Wingspan in feet (ft)
    gamma (float): Taper ratio (ratio of tip chord to root chord)
    t_ct (float): Tip thickness-to-chord ratio
    t_cr (float): Root thickness-to-chord ratio
    
    Returns:
    float: Estimated fuel volume in cubic feet (ft³)
    
    Formula:
    The calculation is based on an empirical estimation:
    V_WF = 0.54 * (S² / b) * (t_cr) * ((1 + gamma * sqrt(tau) + (gamma²) * tau) / ((1 + gamma)²))
    where tau = t_ct / t_cr
    """
    
    tau = t_ct / t_cr  # Ratio of tip thickness-to-chord ratio to root thickness-to-chord ratio
    
    # Compute fuel volume capacity
    V_WF = 0.54 * (S**2 / b) * (t_cr) * ((1 + gamma * (tau**0.5) + (gamma**2) * tau) / ((1 + gamma)**2))
    
    return V_WF


def CalcMAC(c_root, taper):
    """
    Calculates the Mean Aerodynamic Chord (MAC) for a trapezoidal wing.
    
    Parameters:
    c_root (float): Root chord length in feet (ft)
    taper (float): Taper ratio (tip chord / root chord)
    
    Returns:
    float: Mean Aerodynamic Chord (MAC) in feet (ft)
    
    Formula:
    MAC = (2/3) * c_root * ((1 + taper + taper²) / (1 + taper))
    """
    
    MAC = (2/3) * c_root * ((1 + taper + taper**2) / (1 + taper))
    return MAC



# Example usage
S = 193.6 # Wing area in square feet of fore wing in ft^2
b = 44     # Wingspan in feet
gamma = 0.7 # Taper ratio
t_ct = 0.18 # Tip thickness-to-chord ratio
t_cr = 0.18 # Root thickness-to-chord ratio
c_root = 5.18 # root chord length in ft
c_root_rudder = 5.12
gamma_rudder = 0.8

# Compute fuel volume
V_WF = CalcFuelVolume(S, b, gamma, t_ct, t_cr)

# Compute MAC
MAC = CalcMAC(c_root, gamma)

# Display the results
print(f"Estimated wing fuel volume: {V_WF:.2f} ft³")
print(f"Mean Aerodynamic Chord (MAC): {MAC:.2f} ft")

MAC = CalcMAC(c_root_rudder, gamma_rudder)
print(f"Mean Aerodynamic Chord (Rudder): {MAC:.2f} ft")
