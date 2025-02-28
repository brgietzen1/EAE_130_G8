import numpy as np

#distance from CG back to tail
l_h = 
#horizontal tail reference area
S_h = 
#mean aerodynamic chord
chord = 
#wing reference area
S_w = 

#wing lift coefficient
#aspect ratio
AR = 

eta = .97
#wing sweep
delta = 
#mach number
M =

C_lw = 2*np.pi*AR/(2 + np.sqrt((AR/eta)**2 * (1+np.tan(delta**2 - M**2)) + 4))

neutral_point = (l_h * S_h / chord / S_w) * C_lh * (C_lw)^(-1) - C_mf * (C_lw)^(-1) 