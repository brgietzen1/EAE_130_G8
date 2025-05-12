import numpy as np

def quarter_chord_distance(x_w, x_h, c_r_w, c_r_h, taper_w, taper_h, sweep_w, sweep_h, b_w, b_h):
    """
    Computes the quarter-chord positions of a swept wing and a swept horizontal stabilizer,
    and calculates their longitudinal distance.
    
    Parameters:
        x_w: float - Leading-edge root position of the wing.
        x_h: float - Leading-edge root position of the horizontal stabilizer.
        c_r_w: float - Root chord of the wing.
        c_r_h: float - Root chord of the horizontal stabilizer.
        taper_w, taper_h: float - Taper ratio (tip chord/root chord) of the wing and stabilizer.
        sweep_w, sweep_h: float - Leading-edge sweep angle (degrees) of the wing and stabilizer.
        b_w, b_h: float - Span of the wing and stabilizer.
        
    Returns:
        tuple - Quarter-chord positions of the wing and stabilizer, and their longitudinal distance.
    """
    # Convert sweep angles to radians
    sweep_w = np.radians(sweep_w)
    sweep_h = np.radians(sweep_h)
    
    # Compute quarter-chord positions
    x_qc_w = x_w + c_r_w / 4 + ((b_w / 2) * np.tan(sweep_w)) / (3 * (1 + taper_w))
    x_qc_h = x_h + c_r_h / 4 + ((b_h / 2) * np.tan(sweep_h)) / (3 * (1 + taper_h))
    
    # Compute quarter-chord distance
    delta_x_qc = x_qc_h - x_qc_w
    
    return x_qc_w, x_qc_h, delta_x_qc

# Given data
x_w = 9.0   # ft (Root leading-edge of the wing)
x_h = 31.5  # ft (Root leading-edge of the stabilizer)
c_r_w = 5.18
c_r_h = 5.18
taper_w = 0.7
taper_h = 0.7
sweep_w = 15  # degrees
sweep_h = -15 # degrees
b_w = 30  # ft (Example wing span)
b_h = 12  # ft (Example stabilizer span)

# Compute the quarter-chord positions and distance
x_qc_w, x_qc_h, delta_x_qc = quarter_chord_distance(x_w, x_h, c_r_w, c_r_h, taper_w, taper_h, sweep_w, sweep_h, b_w, b_h)

# Print results
print(f"Quarter-chord position (Wing): {x_qc_w:.4f} ft")
print(f"Quarter-chord position (Stabilizer): {x_qc_h:.4f} ft")
print(f"Quarter-chord distance: {delta_x_qc:.4f} ft")