import numpy as np

#distance from CG back to tail
l_h = 17
#horizontal tail reference area
S_h = 162.5
#mean aerodynamic chord
chord = 4
#wing reference area
S_w = 162.5

#wing lift coefficient
#wing aspect ratio
AR_w = 8.5
#constance from slides
eta = .97
#wing sweep
delta = 15.5
#mach number
M = .16

#tail lift coefficient
#tail aspect ratio
AR_h = 8.5
#tail sweep
delta_h = -15.5

#fuselage pitching moment
#constant depending on wing 1/4 chord position
K_f = .115
#fuselage width
w_f = 3.79
#fuselage length
L_f = 35.5

C_l_w = 2*np.pi*AR_w/(2 + np.sqrt((AR_w/eta)**2 * (1+np.tan(delta**2 - M**2)) + 4))

deda = 2/(np.pi*AR_w) * C_l_w

C_l_h_clean = 2*np.pi*AR_h/(2 + np.sqrt((AR_h/eta)**2 * (1+np.tan(delta_h**2 - M**2)) + 4))

C_l_h = C_l_h_clean * (1-deda)

C_m_fus = K_f * w_f**2 * L_f / (S_w * chord) * (C_l_w)**(-1)

neutral_point = (l_h * S_h / chord / S_w) * C_l_h * (C_l_w)**(-1) - C_m_fus * (C_l_w)**(-1)

print(neutral_point)