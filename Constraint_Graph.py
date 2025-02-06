import numpy as np
import matplotlib.pyplot as plt

#wing loading vector
wing_loading = np.linspace(10, 50, 100)
#stall velocity
v_stall = 
#maximum coefficient of lift at cruise
cl_max_cruise = 
#Takeoff distance
s_TO = 
#density on the ground (should be at sea level for Davis)
rho_ground = 
#density at sea level
rho_sea = 
#maximum coefficient of lift at takeoff
cl_max_TO = 
#total landing distance
s_L = 
#obstacle clearance distance
s_a = 
#maximum coefficient of lift at landing conditions
cl_max_land = 
e = 
#Aspect ratio
AR =
k_s = 
cd_0 = 
cl_max_climb = 


#Stall Constraint (plot on x axis as vertical line) !NEEDS POWER CORRECTION
def stall_constraint(v_stall, cl_max_cruise):
    wing_loading_stall = 1/2 * v_stall**2 * cl_max_cruise
    return wing_loading_stall

#Takeoff constraint (function of W/S) !NEEDS POWER CORRECTION
def takeoff_constraint(s_TO, rho_ground, rho_sea, cl_max_TO, wing_loading):
    a = .0149
    b = 8.314
    c = -s_TO
    coefficients = [a, b, c]
    solutions = np.roots(coefficients)
    TOP_23 = max(solutions)
    takeoff_constraint = wing_loading / (rho_ground/rho_sea * cl_max_TO * TOP_23)
    return takeoff_constraint

#Landing constraint (plot on x axis as multiple vertical lines for each cl) !NEEDS POWER CORRECTION
def landing_constraint(s_L, s_a, cl_max_land, rho_ground, rho_sea):
    wing_loading_landing = (s_L - s_a) / 80 * cl_max_land * rho_ground/rho_sea
    return wing_loading_landing

#Climb constraint (plot on y axis as horizontal line) !NEEDS POWER CORRECTION
def climb_constraint(e, AR, k_s, cd_0, cl_max_climb):
    G = 9.8
    k = 1 / (np.pi* e * AR)
    climb_constraint = k_s**2 * cd_0 / cl_max_climb + k * cl_max_climb / k_s**2 + G
    climb_constraint_corrected = 1/.8 * 1/.94 * climb_constraint
    return climb_constraint_corrected

#Cruise constraint
def cruise_constraint(e, AR, W_cruise, W_takeoff, wing_loading, v_cr, eta_p, q_cr, cd_0, P_cruise, P_takeoff):
    k = 1 / (np.pi* e * AR)
    wing_loading_cruise = wing_loading * (W_cruise / W_takeoff)
    cruise_constraint = v_cr / eta_p * (q_cr * cd_0 / (wing_loading_cruise) + k / q_cr * wing_loading_cruise)
    cruise_constraint_corrected = W_cruise / W_takeoff / (P_cruise / P_takeoff) * cruise_constraint
    return cruise_constraint_corrected



