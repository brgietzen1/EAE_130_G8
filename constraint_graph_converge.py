import numpy as np
import Constraint_Graph as cg
S = np.linspace(100, 500, 200)
for i in range(len(S)):
    S_0 = S(i)
    P(i) = P_guess
    tolerance = .01
    converged = False
    while converged == False:
        W = f(S_0, T(i))
        wing_loading = W/S_0
        power_loading_new = f(wing_loading)
        P_new = 1/power_loading_new * W
        if abs(P_new - T(i)) <= tolerance:
            converged = True
        P(i) = P_new

print(P)
    
