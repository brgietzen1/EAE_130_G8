import numpy as np

#distance from aerodynamic center to back tail
l_h = 21.3967
#horizontal tail reference area
S_h = 193.6
#mean aerodynamic chord
chord = 4.45
#wing reference area
S_w = 193.6

#wing lift coefficient
#wing aspect ratio
AR_w = 10
#constance from slides
eta = .97
#wing sweep
delta = np.deg2rad(5)
#mach number
M = 180.446/1116.133

#tail lift coefficient
#tail aspect ratio
AR_h = 10
#tail sweep
delta_h = np.deg2rad(-5)

#fuselage pitching moment
#constant depending on wing 1/4 chord position
K_f = .344
#fuselage width
w_f = 4.4
#fuselage length
L_f = 35.5

C_l_w = 2*np.pi*AR_w/(2 + np.sqrt((AR_w/eta)**2 * (1+np.tan(delta**2 - M**2)) + 4))

deda = 2/(np.pi*AR_w) * C_l_w 

C_l_h_clean = 2*np.pi*AR_h/(2 + np.sqrt((AR_h/eta)**2 * (1+np.tan(delta_h**2 - M**2)) + 4))

C_l_h = C_l_h_clean * (1-deda) 

C_m_fus = K_f * w_f**2 * L_f / (S_w * chord) * (C_l_w)**(-1)

neutral_point = (l_h * S_h / chord / S_w) * C_l_h * (C_l_w)**(-1) - C_m_fus * (C_l_w)**(-1)

print(neutral_point)
print(neutral_point + 9 + 5.18/4)
print( 9 + 5.18/4 + neutral_point - .4*chord)