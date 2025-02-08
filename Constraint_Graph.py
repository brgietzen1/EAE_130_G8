import numpy as np
import matplotlib.pyplot as plt

#!!! note that all power loading values are in lbf/hp
#wing loading vector [lbf/ft^2]
wing_loading = np.linspace(10, 50, 100)
#stall velocity [ft/s]
v_stall = 133
#overall maximum coefficient of lift
cl_max = 2.5
#maximum coefficient of lift at cruise
cl_max_cruise = 1.2
#Takeoff distance [ft]
s_TO = 800
#density on the ground (should be at sea level for Davis) [slug/ft^3]
rho_ground = .002378
#density at sea level [slug/ft^3]
rho_sea = .002378
#density at cruise [slug/ft^3]
rho_cruise = .002378
#maximum coefficient of lift at takeoff
cl_max_TO = 2.14
#total landing distance [ft]
s_L = 1640
#obstacle clearance distance [ft]
s_a = 450
#maximum coefficient of lift at landing conditions
cl_max_land = 2.5
#span efficiency factor
e = 0.85
#Aspect ratio
AR = 9.22
#stall factor (in FAR 23)
k_s = 1.3
#drag coefficient at zero lift
cd_0 = .02862
#max coefficient of lift at climb
cl_max_climb = 1.8
#cruise weight [lbm]
W_cruise = 12000
#takeoff weight [lbm]
W_takeoff = 12214
#cruise speed [ft/s]
v_cruise = 211
#propeller efficiency
eta_p = 0.8
#dynamic pressure at cruise [lbf/ft^2] slugs/ft^3 * ft^2/s^2 = slugs/ft/s^2 = lbf/s^2/ft^2
q_cruise = (1/2*rho_cruise * v_cruise**2)
#power at cruise [hp]
P_cruise = 500
#power at takeoff [hp]
P_takeoff = 600
#turn radius [ft]
R_turn = 400


#Stall Constraint (plot on x axis as vertical line)
def stall_constraint(v_stall, rho_ground, cl_max):
    wing_loading_stall = 1/2 * rho_ground * v_stall**2 * cl_max
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
    return takeoff_constraint

#Landing constraint (plot on x axis as multiple vertical lines for each cl)
def landing_constraint(s_L, s_a, cl_max_land, rho_ground, rho_sea):
    wing_loading_landing = (s_L - s_a) / 80 * cl_max_land * rho_ground/rho_sea
    return wing_loading_landing

#Climb constraint (plot on y axis as horizontal line) !POWER CORRECTION
def climb_constraint(e, AR, k_s, cd_0, cl_max_climb, rho_ground, eta_p):
    G = 0.012
    k = 1 / (np.pi* e * AR)
    climb_constraint = k_s**2 * cd_0 / cl_max_climb + k * cl_max_climb / k_s**2 + G
    climb_constraint_corrected = 1/.8 * 1/.94 * climb_constraint
    V = np.sqrt(1/rho_ground / cl_max_climb * wing_loading)
    climb_constraint_corrected_power = 550*eta_p/V * 1/climb_constraint_corrected
    return climb_constraint_corrected

#Cruise constraint (function of W/S) !POWER CORRECTION
def cruise_constraint(e, AR, W_cruise, W_takeoff, wing_loading, v_cruise, eta_p, q_cruise, cd_0, P_cruise, P_takeoff, rho_cruise, cl_max_cruise):
    k = 1 / (np.pi* e * AR)
    wing_loading_cruise = wing_loading * (W_cruise / W_takeoff)
    cruise_constraint = v_cruise / 550 /  eta_p * (q_cruise * cd_0 / (wing_loading_cruise) + k / q_cruise * wing_loading_cruise)
    cruise_constraint_corrected = 1/(W_cruise / W_takeoff / (P_cruise / P_takeoff) * cruise_constraint)
    return cruise_constraint_corrected

#Ceiling constraint (plot on y axis as horizontal line) !POWER CORRECTION
def ceiling_constraint(e, AR, cd_0, rho_cruise, cl_max_cruise, wing_loading, eta_p):
    k = 1 / (np.pi* e * AR)
    ceiling_constraint = 2*np.sqrt(k*cd_0)
    V = np.sqrt(2/rho_cruise/cl_max_cruise * wing_loading)
    ceiling_constraint_power = 550*eta_p/V * 1/ceiling_constraint
    return ceiling_constraint

#Manuever constraint (function of W/S) !POWER CORRECTION
def manuever_constraint(e, AR, v_cruise, R_turn, q_cruise, cd_0, wing_loading, rho_ground, cl_max_climb, eta_p):
    k = 1 / (np.pi* e * AR)
    g = 32.2
    n = np.sqrt((v_cruise**2 / R_turn / g)**2 + 1)
    maneuver_constraint = q_cruise*cd_0/wing_loading + k * n**2 / q_cruise * wing_loading
    V = np.sqrt(2/rho_ground / cl_max_climb * wing_loading)
    maneuver_constraint_power = 550*eta_p/V * 1/maneuver_constraint
    return maneuver_constraint

stall_constraint = stall_constraint(v_stall, rho_ground, cl_max)
print(stall_constraint)
'''plt.figure()
plt.plot(wing_loading, wing_loading)
plt.axvline(x=stall_constraint, color='r', linestyle='--', label="x = 2")
plt.title("stall")
plt.show()'''


takeoff_constraint = takeoff_constraint(s_TO, rho_ground, rho_sea, cl_max_TO, wing_loading, eta_p)
print(takeoff_constraint)
'''plt.figure()
plt.plot(wing_loading, takeoff_constraint)
plt.title("takeoff")
plt.show()'''

landing_constraint = landing_constraint(s_L, s_a, cl_max_land, rho_ground, rho_sea)
print(landing_constraint)
'''plt.figure()
plt.plot(wing_loading, takeoff_constraint)
plt.axvline(x=landing_constraint, color='g', linestyle='--', label="x = 2")
plt.title("landing")
plt.show()'''

climb_constraint = climb_constraint(e, AR, k_s, cd_0, cl_max_climb, rho_ground, eta_p)
print(climb_constraint)
'''plt.figure()
plt.plot(wing_loading, climb_constraint)
plt.title("climb")
plt.show()'''

'''cruise_constraint = cruise_constraint(e, AR, W_cruise, W_takeoff, wing_loading, v_cruise, eta_p, q_cruise, cd_0, P_cruise, P_takeoff, rho_cruise, cl_max_cruise)
print(cruise_constraint)'''
'''plt.figure()
plt.plot(wing_loading, cruise_constraint)
plt.title("cruise")
plt.show()'''

ceiling_constraint = ceiling_constraint(e, AR, cd_0, rho_cruise, cl_max_cruise, wing_loading, eta_p)
print(ceiling_constraint)
'''plt.figure()
plt.plot(wing_loading, ceiling_constraint)
plt.title("ceiling")
plt.show()'''

manuever_constraint = manuever_constraint(e, AR, v_cruise, R_turn, q_cruise, cd_0, wing_loading, rho_ground, cl_max_climb, eta_p)
print(manuever_constraint)
'''plt.figure()
plt.plot(wing_loading, manuever_constraint)
plt.title("maneuver")
plt.show()'''

'''plt.figure()
plt.plot(wing_loading, takeoff_constraint, label='takeoff')
plt.plot(wing_loading, climb_constraint, label='climb')
plt.plot(wing_loading, cruise_constraint, label='cruise')
plt.plot(wing_loading, ceiling_constraint, label='ceiling')
plt.plot(wing_loading, manuever_constraint, label='maneuver')
plt.axvline(x=stall_constraint, color='r', linestyle='--', label='stall')
plt.axvline(x=landing_constraint, color='g', linestyle='--', label='landing')
plt.legend()
plt.title("All")
plt.xlabel("W/S [lbf/ft^2]")
plt.ylabel("W/P [lbf/hp]")
plt.show()'''

plt.figure()
plt.plot(wing_loading, takeoff_constraint, label='takeoff')
plt.plot(wing_loading, np.ones_like(wing_loading)*climb_constraint, label='climb')
'''plt.plot(wing_loading, cruise_constraint, label='cruise')'''
plt.plot(wing_loading, np.ones_like(wing_loading)*ceiling_constraint, label='ceiling')
plt.plot(wing_loading, manuever_constraint, label='maneuver')
plt.axvline(x=stall_constraint, color='r', linestyle='--', label='stall')
plt.axvline(x=landing_constraint, color='g', linestyle='--', label='landing')
plt.legend()
plt.title("All")
plt.xlabel("W/S [lbf/ft^2]")
plt.ylabel("W/P [lbf/hp]")
plt.show()