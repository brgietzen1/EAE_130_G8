import numpy as np

# calculator helps find the x_coordinate of the leading edge of a airfoil section in a tapered wing

def chord_at_y(y, semi_span, c_root, c_tip):
    return c_root - ((c_root - c_tip) / semi_span) * y

def x_le_from_root_le(y, semi_span, sweep_qc_deg, c_root, c_tip):
    gamma_rad = np.radians(sweep_qc_deg)
    c_y = chord_at_y(y, semi_span, c_root, c_tip)
    x_le = c_root / 4 + y * np.tan(gamma_rad) - c_y / 4
    return x_le


semi_span = 8.9  # ft
sweep_qc_deg = 15.0  # degrees
c_root = 5.00  # ft
c_tip = 4.00  # ft

y = [0, 1.37, 7.5, 8.9]  # ft from root

for y_val in y:
    x_le = x_le_from_root_le(y_val, semi_span, sweep_qc_deg, c_root, c_tip)
    print(f"Leading edge x at y = {y} ft: {x_le:.4f} ft")

for y_val in y:

    c_y = chord_at_y(y_val, semi_span, c_root, c_tip)
    print(f"Chord lenght at y = {y} ft: {c_y:.4f} ft")