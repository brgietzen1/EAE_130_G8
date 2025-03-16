import numpy as np


W_payload = 4000
W_pax = 190
#takeoff weight
W_TO =
#wing aspect ratio
A = 10
#wing quarter-chord sweep angle
del_quar =
#design ultimate load factor
n_ult =
#wing area [ft^2]
S =
#wing taper ratio
lam =
#maximum wing thickness ratio
tcm =
#maximum level speed at sealevel [kts]
V_H =

#vertical tail area [ft]
S_v = 
#vertical tail span [ft]
b_v =
#vertical tail maximum root thickness [ft]
trv =

#fuselage length [ft]
l_f = 35.5
#max fuselage width [ft]
w_f = 4.4
#max fuselage height [ft]
h_f = 6.5
#design cruise speed [KEAS]
V_C =

#shock strut length for main gear [ft]
l_s_m = 
#design landing weight [lbs]
W_L = W_TO-W_payload

#engine takeoff power [shp]
P_T = 620
#number of engine
N_e = 1

#mission fuel weight [lbs]
W_F = 
#constant
K_fsp = 5.87
#fuel fraction of tanks which are integral
fracfuel = 
#number of separate engine tanks
N_t = 

#number of passengers 
N_pax = 1

#design dive mach number
M_D =





def wing_weight(W_TO, n_ult, A, del_quar, S, lam, tcm, V_H):
    wing_weight = 96.948*(W_TO*n_ult/(10**5)**(.65)*(A/np.cos(del_quar))**(.57)*(S/100)**(.61)**((1+lam)/2*(tcm))**(.36)*(1+V_H/500)**(.5))**(.993)
    return wing_weight

def empennage_weight(W_TO, n_ult, S_v, b_v, trv):
    W_V = 98.5*((W_TO*n_ult/(10**5))**(.87)*(S_v/100)**(1.2)*.289*(b_v/trv)**(.5))**(.458)
    empennage_weight = W_V
    return empennage_weight

def fuselage_weight(W_TO, n_ult, l_f, w_f, h_f, V_C):
    fuselage_weight = 200*((W_TO*n_ult/(10**5))**(.286)*(l_f/10)**(.857)*((w_f+h_f)/10)*(V_C/100)**(.338))**(1.1)
    return fuselage_weight

# def nacelle_weight():


def gear_weight(l_s_m, W_L, n_ult_l=5.7):
    gear_weight = 0.054*(l_s_m)**(.501)*(W_L*n_ult_l)**(.684)
    return gear_weight


def engine_weight(P_TO, N_e, k_p=.45):
    #includes air induction, propeller, propulsion system
    W_eng = k_p*P_TO
    engine_weight = 2.575*W_eng**(.992)*N_e
    return engine_weight

def fuel_system_weight(W_F, K_fsp, fracfuel, N_t, N_e):
    W_fs = 2.49*((W_F/K_fsp)**(.6)*(1/(1+fracfuel))**(.3)*(N_t)**(.2)*N_e**(.13))**(1.21)
    return W_fs

def flight_control_sys_weight(W_TO):
    W_fc = 1.08*(W_TO)**(.7)
    return W_fc

def electrical_system_weight(W_fs, W_aie):
    W_els = 426*((W_fs + W_aie)/1000)**(.51)
    return W_els

def avionica_weight(N_pax):
    W_iae = 33*N_pax
    return W_iae

def AC_weight(W_TO, N_pax, W_iae, M_D):
    W_api = .265*(W_TO)**(.52)*(N_pax)**(.68)*(W_iae)**(.17)*M_D**(.08)
    return W_api
    
    
    
# def oxygen_sys_weight

def APU_weight(W_TO):
    W_APU = .0085*W_TO
    return W_APU

def furnishings_weight(N_pax, W_TO):
    W_fur = .412*N_pax**(1.145)*W_TO**(.489)
    return W_fur
