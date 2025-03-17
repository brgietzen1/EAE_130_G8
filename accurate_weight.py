import numpy as np

def wing_weight(W_TO, n_ult, A, del_quar, S, lam, tcm, V_H):
    wing_weight_USAF = 96.948*((W_TO*n_ult/(10**5))**(.65)*(A/np.cos(del_quar))**(.57)*(S/100)**(.61)*((1+lam)/2*(tcm))**(.36)*(1+V_H/500)**(.5))**(.993)
    wing_weight_ces = .04674*(W_TO)**(.397)*(S**.36)*n_ult**.397*A**1.712
    return (wing_weight_USAF + wing_weight_ces) / 2

def empennage_weight(W_TO, n_ult, S_v, b_v, trv, A_v, del_quar_v):
    W_V_USAF = 98.5*((W_TO*n_ult/(10**5))**(.87)*(S_v/100)**(1.2)*.289*(b_v/trv)**(.5))**(.458)
    W_V_tor = .04*(n_ult*S_v**2)**.75
    empennage_weight = (W_V_USAF + W_V_tor )/2
    return empennage_weight

def fuselage_weight(W_TO, n_ult, l_f, w_f, h_f, V_C, N_pax):
    fuselage_weight_USAF = 200*((W_TO*n_ult/(10**5))**(.286)*(l_f/10)**(.857)*((w_f+h_f)/10)*(V_C/100)**(.338))**(1.1)
    fuselage_weight_ces = .04682*W_TO**.692*N_pax*.374*l_f**.59/100
    return (fuselage_weight_USAF + fuselage_weight_ces) / 2


def engine_weight(P_TO, N_e, N_p, N_bl, D_p, k_p=.55):
    #includes nacelle, air induction, propeller, propulsion system
    W_eng = k_p*P_TO
    engine_weight_USAF = 2.575*W_eng**(.992)*N_e
    W_n = .37*P_TO 
    W_prop = 31.92*N_p*N_bl**.391*((D_p)*P_TO/N_e/1000)**.782
    W_p = 65
    engine_weight_ces = W_eng + W_n + W_prop + W_p
    return (engine_weight_USAF + engine_weight_ces)/2

def fuel_system_weight(W_F, K_fsp, fracfuel, N_t, N_e):
    W_fs = 2.49*((W_F/K_fsp)**(.6)*(1/(1+fracfuel))**(.3)*(N_t)**(.2)*N_e**(.13))**(1.21)
    return W_fs

def flight_control_sys_weight(W_TO):
    W_fc = 1.08*(W_TO)**(.7)
    return W_fc

def electrical_system_weight(W_fs, W_iae):
    W_els = 426*((W_fs + W_iae)/1000)**(.51)
    return W_els

def avionics_weight(N_pax):
    W_iae = 33*N_pax
    return W_iae

def AC_weight(W_TO, N_pax, W_iae, M_D):
    W_api = .265*(W_TO)**(.52)*(N_pax)**(.68)*(W_iae)**(.17)*M_D**(.08)
    return W_api

def APU_weight(W_TO):
    W_APU = .013*W_TO
    return W_APU

def furnishings_weight(N_pax, W_TO):
    W_fur = .412*N_pax**(1.145)*W_TO**(.489)
    return W_fur


def total_weight_calc():
    W_TO = 10910
    #payload weight
    W_payload = 4000
    #passenger weight
    W_pax = 190
    #landing gear weight
    W_lg = 384 + 53
    #boom weight
    W_boom = 150
    #mission fuel weight [lbs]
    W_F = 642.34
    #wing aspect ratio
    A = 10
    #wing quarter-chord sweep angle
    del_quar = np.deg2rad(5)
    #wing area [ft^2]
    S = 193.6
    #wing taper ratio
    lam = 0.7
    #maximum wing thickness ratio
    tcm = 0.18
    #maximum level speed at sealevel [kts]
    V_H = 394/309*117
    #vertical tail area [ft]
    S_v = 40.4624
    #vertical tail span [ft]
    b_v = 8.9
    #vertical tail maximum root thickness [ft]
    trv = 0.16
    #vertical tail aspect ratio
    A_v = 1.88
    #vertical tail quarter chord sweep angle
    del_quar_v = np.deg2rad(15)
    #fuselage length [ft]
    l_f = 35.5
    #max fuselage width [ft]
    w_f = 4.4
    #max fuselage height [ft]
    h_f = 6.5
    #design cruise speed [KEAS]
    V_C = 128
    #engine takeoff power [shp]
    P_TO = 620
    #number of engine
    N_e = 1
    #number of propeller
    N_p = 1
    #number of propeller blades
    N_bl = 5
    #diameter of propeller
    D_p = 9
    #constant
    K_fsp = 5.87
    #fuel fraction of tanks which are integral
    fracfuel = 1
    #number of separate engine tanks
    N_t = 1
    #number of passengers 
    N_pax = 1
    #design dive mach number
    M_D = .438
    tolerance = .05
    iter = 0
    W_TO_new = 0
    W_TO_old = W_TO
    while abs(W_TO_new-W_TO_old)/W_TO_old >= .05:
        W_TO_old = W_TO
        #positive design limit load factor
        n_lim_pos = 2.1 + 24000/(W_TO + 10000)
        #design ultimate load factor
        n_ult = 1.5*n_lim_pos
        #design landing weight [lbs]
        W_L = W_TO-W_payload
        W_w = wing_weight(W_TO, n_ult, A, del_quar, S, lam, tcm, V_H)
        W_em = empennage_weight(W_TO, n_ult, S_v, b_v, trv, A_v, del_quar_v)
        W_fus = fuselage_weight(W_TO, n_ult, l_f, w_f, h_f, V_C, N_pax)
        W_eng = engine_weight(P_TO, N_e, N_p, N_bl, D_p, k_p=.55)
        W_fs = fuel_system_weight(W_F, K_fsp, fracfuel, N_t, N_e)
        W_fc = flight_control_sys_weight(W_TO)
        W_iae = avionics_weight(N_pax)
        W_els = electrical_system_weight(W_fs, W_iae)
        W_api = AC_weight(W_TO, N_pax, W_iae, M_D)
        W_APU = APU_weight(W_TO)
        W_fur = furnishings_weight(N_pax, W_TO)
        W_winglets = 1/12.7*2*W_w

        W_TO_new = W_payload + W_F + W_pax + W_boom + 2*W_w + W_em + W_fus + W_lg + W_eng + W_fs + W_fc + W_iae + W_els + W_api + W_APU + W_fur + W_winglets
        W_TO = W_TO_new
        iter = iter + 1

    return iter, W_payload, W_F, W_pax, n_lim_pos, W_w, W_em, W_fus, W_lg, W_eng, W_fs, W_fc, W_iae, W_els, W_api, W_APU, W_fur, W_winglets, W_TO_new
    
    
iter, W_payload, W_F, W_pax, n_lim_pos, W_w, W_em, W_fus, W_lg, W_eng, W_fs, W_fc, W_iae, W_els, W_api, W_APU, W_fur, W_winglets, W_TO_new = total_weight_calc()
print("Iterations: ", iter)
print("Payload: ", W_payload)
print("Fuel: ", W_F)
print("Passenger: ", W_pax)
print("Positive Load Limit:", n_lim_pos)
print("Wings: ", 2*W_w)
print("Winglets: ", W_winglets)
print("Empennage: ", W_em)
print("Fuselage: ", W_fus)
print("Landing Gear: ", W_lg)
print("Engine: ", W_eng)
print("Fuel Systems: ", W_fs)
print("Flight Control Systems: ", W_fc)
print("Avionics: ", W_iae)
print("Electrical Systems: ", W_els)
print("AC, Pressure, Anti-Ice Systems: ", W_api)
print("APU: ", W_APU)
print("Furnishings: ", W_fur)
print("Gross Takeoff Weight: ", W_TO_new)


