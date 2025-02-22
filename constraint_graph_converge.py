import numpy as np
import matplotlib.pyplot as plt
import Constraint_Graph as cg
import a4_weight_estimation as we

S = np.linspace(100, 500, 400)
P = np.zeros_like(S)
for i in range(len(S)):
    S_0 = S[i]
    P[i] = 300
    tolerance = .01
    converged = False
    while converged == False:
        results = we.solve_takeoff_weight_2(S_0, P[i], S_design = 12214.12 / 31.25, P_design = 12214.12 / 22.2)
        takeoff_weight, empty_weight_frac, fuel_weight_frac,\
        empty_weight, fuel_weight, iterations = results
        W = takeoff_weight
        wing_loading = W/S_0
        takeoff_constraint, takeoff_constraint_power = cg.takeoff_constraint(s_TO=1000, rho_ground=.002378, rho_sea=.002378, cl_max_TO=2.23, wing_loading=wing_loading, eta_p=0.9)
        power_loading_new = takeoff_constraint_power
        P_new = 1/power_loading_new * W
        if abs(P_new - P[i]) <= tolerance:
            converged = True
        P[i] = P_new

print(S)
print(P)
    
plt.figure()
plt.plot(S, P, label='Takeoff')
plt.legend()
plt.grid()
plt.title("Power (P) vs. Wing Area (S)")
plt.xlabel("S [ft\u00b2]")
plt.ylabel("P [hp]")
plt.show()