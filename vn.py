import numpy as np
import matplotlib.pyplot as plt

V = np.linspace(0, 400, 400)
def max_load_limit_factor(W):
    n_max = 2.1 + 24000/(W + 10000)
    return n_max

def neg_load_limit_factor(n_max):
    n_neg = -0.4*n_max
    return n_neg

def stall_load(rho, rho_SL, V, CL_max, W, S):
    n_stall = rho_SL*(V)**2*CL_max/2/(W/S)
    return n_stall

def stall_speed(W, rho, rho_SL, S, CL_max):
    v_s1 = np.sqrt(2*W/rho/S/CL_max)*np.sqrt(rho/rho_SL)
    return v_s1

def stall_speed_neg(W, rho, rho_SL, S, CL_min):
    v_sneg1 = np.sqrt(-2*W/rho/S/CL_min)*np.sqrt(rho/rho_SL)
    return v_sneg1 

def corner_speed(n_max, W, S, rho, rho_SL, CL_max):
    V_A = np.sqrt(n_max*2*W/S/rho_SL/CL_max)
    return V_A

def max_level_flight_speed(V_D, rho, rho_SL):
    V_H = 0.9*V_D
    return V_H

def dive_speed(V_C, rho, rho_SL):
    V_D = 1.25*V_C
    return V_D

W = 8998
rho = 1.758e-3
rho_SL = 2.38e-3
S = 387
CL_max = 1.612
CL_min = -.49

n_max = max_load_limit_factor(W)
n_neg = neg_load_limit_factor(n_max)
n_stall = stall_load(rho, rho_SL, V, CL_max, W, S)
n_stall_neg = stall_load(rho, rho_SL, V, CL_min, W, S)

v_s1 = stall_speed(W, rho, rho_SL, S, CL_max)
v_sneg1 = stall_speed_neg(W, rho, rho_SL, S, CL_min)
V_C = 196*np.sqrt(rho/rho_SL)
V_A = corner_speed(n_max, W, S, rho, rho_SL, CL_max)
V_D = dive_speed(V_C, rho, rho_SL)
V_H = max_level_flight_speed(V_D, rho, rho_SL)

plt.figure()
plt.plot(V, n_stall, color='blue')
plt.plot(V, n_stall_neg, color='blue')
plt.axvline(v_s1, ymin=0.429, ymax=.571, color='red', label='V_s1')
plt.axvline(v_sneg1, ymin=0.286, ymax = .429, color='skyblue', label='V_s-1')
plt.axvline(V_C, ymin=.33, ymax=.76, color='green', label='V_C')
plt.axvline(V_A, ymin=0.429,ymax=(3+n_max)/7, color='cyan', label='V_A')
plt.axvline(V_D, ymin=0.429, ymax=(n_max+3)/7, color='magenta', label='V_D')
plt.axvline(V_H, ymin=.3, ymax=.85, color='yellow', label='V_H')
plt.hlines(n_max, xmin=V_A, xmax=400, color='blue')
plt.hlines(n_neg, xmin=233, xmax=400, color='blue')
plt.hlines(0, xmin=0, xmax=250, color='black')
#plt.plot(V, (0-n_neg)/(V_D-V_H)*(V-V_D), color='blue')

plt.legend()


plt.axis([0, 250, -3, 4])



plt.show()





