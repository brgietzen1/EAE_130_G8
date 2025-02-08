import numpy as np
import matplotlib.pyplot as plt

#!!! note that all power loading values are in lbf/hp
#wing loading vector [lbf/ft^2]
wing_loading = np.linspace(10, 50, 100)
#stall velocity [ft/s]
v_stall = 
#overall maximum coefficient of lift
cl_max =
#maximum coefficient of lift at cruise
cl_max_cruise = 
#Takeoff distance [ft]
s_TO = 
#density on the ground (should be at sea level for Davis) [lbm/ft^3]
rho_ground = 
#density at sea level [lbm/ft^3]
rho_sea = 
#density at cruise [lbm/ft^3]
rho_cruise = 
#maximum coefficient of lift at takeoff
cl_max_TO = 
#total landing distance [ft]
s_L = 
#obstacle clearance distance [ft]
s_a = 
#maximum coefficient of lift at landing conditions
cl_max_land = 
#span efficiency factor
e = 
#Aspect ratio
AR =
#stall factor (in FAR 23)
k_s = 
#drag coefficient at zero lift
cd_0 = 
#max coefficient of lift at climb
cl_max_climb = 
#cruise weight [lbm]
W_cruise =
#takeoff weight [lbm]
W_takeoff = 
#cruise speed [ft/s]
v_cruise = 
#propeller efficiency
eta_p = 
#dynamic pressure at cruise [lbf/ft^2] lbm/ft^3 * ft^2/s^2 = lmb/ft/s^2/32.2 = lbf/s^2/ft^2
q_cruise = (1/2*rho_cruise * v_cruise**2)/32.2
#power at cruise [hp]
P_cruise =
#power at takeoff [hp]
P_takeoff = 
#turn radius [ft]
R_turn = 


#Stall Constraint (plot on x axis as vertical line)
def stall_constraint(v_stall, cl_max):
    wing_loading_stall = 1/2 * v_stall**2 * cl_max
    return wing_loading_stall

#Takeoff constraint (function of W/S) !POWER CORRECTION
def takeoff_constraint(s_TO, rho_ground, rho_sea, cl_max_TO, wing_loading, eta_p):
    a = .0149
    b = 8.314
    c = -s_TO
    coefficients = [a, b, c]
    solutions = np.roots(coefficients)
    TOP_23 = max(solutions)
    takeoff_constraint = wing_loading / (rho_ground/rho_sea * cl_max_TO * TOP_23)
    V = np.sqrt(2 / rho_ground / cl_max_TO * wing_loading)
    takeoff_constraint_power = 550*eta_p/V * 1/takeoff_constraint
    return takeoff_constraint_power

#Landing constraint (plot on x axis as multiple vertical lines for each cl)
def landing_constraint(s_L, s_a, cl_max_land, rho_ground, rho_sea):
    wing_loading_landing = (s_L - s_a) / 80 * cl_max_land * rho_ground/rho_sea
    return wing_loading_landing

#Climb constraint (plot on y axis as horizontal line) !POWER CORRECTION
def climb_constraint(e, AR, k_s, cd_0, cl_max_climb, rho_ground, eta_p):
    G = 32.2
    k = 1 / (np.pi* e * AR)
    climb_constraint = k_s**2 * cd_0 / cl_max_climb + k * cl_max_climb / k_s**2 + G
    climb_constraint_corrected = 1/.8 * 1/.94 * climb_constraint
    V = np.sqrt(1/rho_ground / cl_max_climb * wing_loading)
    climb_constraint_corrected_power = 550*eta_p/V * 1/climb_constraint_corrected
    return climb_constraint_corrected_power

#Cruise constraint (function of W/S) !POWER CORRECTION
def cruise_constraint(e, AR, W_cruise, W_takeoff, wing_loading, v_cruise, eta_p, q_cruise, cd_0, P_cruise, P_takeoff, rho_cruise, cl_max_cruise):
    k = 1 / (np.pi* e * AR)
    wing_loading_cruise = wing_loading * (W_cruise / W_takeoff)
    cruise_constraint = v_cruise / eta_p * (q_cruise * cd_0 / (wing_loading_cruise) + k / q_cruise * wing_loading_cruise)
    cruise_constraint_corrected = W_cruise / W_takeoff / (P_cruise / P_takeoff) * cruise_constraint
    V = np.sqrt(1/rho_cruise / cl_max_cruise * wing_loading)
    cruise_constraint_corrected_power = 550*eta_p/V * 1/cruise_constraint_corrected
    return cruise_constraint_corrected_power

#Ceiling constraint (plot on vertical axis as horizontal line) !POWER CORRECTION
def ceiling_constraint(e, AR, cd_0, rho_cruise, cl_max_cruise, wing_loading, eta_p):
    k = 1 / (np.pi* e * AR)
    ceiling_constraint = 2*np.sqrt(k*cd_0)
    V = np.sqrt(2/rho_cruise/cl_max_cruise * wing_loading)
    ceiling_constraint_power = 550*eta_p/V * 1/ceiling_constraint
    return ceiling_constraint_power

#Manuever constraint (function of W/S) !POWER CORRECTION
def manuever_constraint(e, AR, v_cruise, R_turn, q, cd_0, wing_loading, rho_ground, cl_max_climb, eta_p):
    k = 1 / (np.pi* e * AR)
    g = 32.2
    n = np.sqrt((v_cruise**2 / R_turn / g)**2 + 1)
    maneuver_constraint = q*cd_0/wing_loading + k * n**2 / q * wing_loading
    V = np.sqrt(2/rho_ground / cl_max_climb * wing_loading)
    maneuver_constraint_power = 550*eta_p/V * 1/maneuver_constraint
    return maneuver_constraint_power

