import numpy as np


def RangeCalculator(eta, c, L_D, W_o, W_f):
    """Simple Breguet Range Equation"""
    R = (eta/c) * (L_D) * np.log(W_o/W_f)

    return R

def SfcUnitConverter(sfc):
    """Converts a specific fuel consumption from lbm/hp/hr to ft^-1
        bsfc: Brake Specific fuel consumption in lbf/hp/hr"""
    sfc_convert = sfc/(1980000)

    return sfc_convert

def LiftCoefficient(Cdo, AR, e):
    """Calculates the lift coefficient when the aircraft is in its clean configuration"""
    k = (np.pi * e* AR)**-1
    C_L = np.sqrt(Cdo/k)
    L_D = 0.94 * C_L / (Cdo + k*(C_L**2))

    return L_D


def PlotPayloadRange(Fuel, Max_Payload, W_e, CD_o, R_ferry, e, eta, AR, c):

    c = SfcUnitConverter(c)

    Fuel_Max_Payload = Fuel * 0.94 # lbf
    Fuel_Reserves = Fuel * 0.06 # lbf
    GTOW = W_e + Fuel_Reserves + Max_Payload + Fuel_Max_Payload

    L_D = LiftCoefficient(CD_o, AR, e)

    # Point A
    R_a = 0
    W_a = Max_Payload

    # Point B
    Wb_o = GTOW
    Wb_f = GTOW - Fuel_Max_Payload
    R_b = RangeCalculator(eta, c, L_D, Wb_o, Wb_f)
    W_b = GTOW

    # Calculating weight of max fuel based on RFP ferry range
    Wd_f = W_e + Fuel_Reserves
    Wd_o = Wd_f * np.exp(R_ferry * c / (eta * L_D))
    max_fuel = Wd_o - Wd_f

    # Point D
    R_d = R_ferry
    W_d =  Wd_o

    # Point C
    Max_Fuel_Payload = GTOW - W_e - Fuel_Reserves - max_fuel
    Wc_o = W_e + Fuel_Reserves + max_fuel + Max_Fuel_Payload
    Wc_f = Wc_o - max_fuel
    R_c = RangeCalculator(eta, c, L_D, Wc_o, Wc_f)
    W_c = Max_Fuel_Payload




# Settings
Total_Mission_Fuel = 822.16
Max_Payload = 3000 # lbf
W_e = 4690 # lbf
CD_o = 0.04048
R_ferry = 3.646e+6 # 600 nmi
e = 1.3567
eta = 0.72 # propeller efficiency
AR = 5
c = 0.6 #lbm/hp/hr

PlotPayloadRange(Total_Mission_Fuel, Max_Payload, W_e, CD_o, R_ferry, e, eta, AR, c)














