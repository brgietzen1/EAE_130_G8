import numpy as np
import matplotlib.pyplot as plt


def RangeCalculator(eta, c, L_D, W_o, W_f):
    """Simple Breguet Range Equation"""
    R = (eta / c) * L_D * np.log(W_o / W_f)  # in ft
    return R

def SfcUnitConverter(sfc):
    """Converts a specific fuel consumption from lbm/hp/hr to ft^-1"""
    return sfc / 1980000

def LiftCoefficient(Cdo, AR, e):
    """Estimates L/D ratio assuming min drag configuration"""
    k = 1 / (np.pi * e * AR)
    C_L = np.sqrt(Cdo / k)
    L_D = 0.94 * C_L / (Cdo + k * C_L ** 2)
    return L_D

def PlotPayloadRange(Fuel, Max_Payload, W_e, CD_o, R_ferry, e, eta, AR, c):
    c = SfcUnitConverter(c)  # convert to ft^-1
    L_D = LiftCoefficient(CD_o, AR, e)
    

    Fuel_Max_Payload = Fuel * 0.94
    Fuel_Reserves = Fuel * 0.06
    GTOW = W_e + Fuel_Reserves + Max_Payload + Fuel_Max_Payload

    # Point A: Zero range, max payload
    R_a = 0
    W_a = Max_Payload

    # Point B: Max payload range 
    Wb_o = GTOW
    Wb_f = GTOW - Fuel_Max_Payload
    R_b = RangeCalculator(eta, c, L_D, Wb_o, Wb_f)
    W_b = Max_Payload

    # Point D: Max range, zero payload (ferry flight)
    Wd_f = W_e + Fuel_Reserves
    Wd_o = Wd_f * np.exp(R_ferry * c / (eta * L_D))
    print(Wd_o)
    max_fuel = Wd_o - Wd_f
    R_d = R_ferry
    W_d = 0

    # Point C: Max range 
    Max_Fuel_Payload = GTOW - W_e - Fuel_Reserves - max_fuel
    Wc_o = W_e + Fuel_Reserves + max_fuel + Max_Fuel_Payload
    Wc_f = Wc_o - max_fuel
    R_c = RangeCalculator(eta, c, L_D, Wc_o, Wc_f)
    W_c = Max_Fuel_Payload

    # Convert range to nautical miles
    ft_to_nmi = 1 / 6076.12
    R_vals = np.array([R_a, R_b, R_c, R_d]) * ft_to_nmi
    W_vals = np.array([W_a, W_b, W_c, W_d])


    print(f"\nMax Payload (lbf): {W_vals[1]:.3f}")
    print(f"Max Payload with Max Fuel (lbf): {W_vals[2]:.3f}")
    print(f"\nMax Fuel For Full Payload (lbf): {Fuel_Max_Payload:.3f}")
    print(f"Max Fuel (lbf): {max_fuel:.3f}")
    print(f"\nMax Payload Range (nmi): {R_vals[1]:.3f}")
    print(f"Max Range With Max Fuel @ MTOW (nmi): {R_vals[2]:.3f}")
    print(f"Max Ferry Range (nmi): {R_vals[3]:.3f}")



    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(R_vals, W_vals, marker='o', linestyle='-', color='b')
    plt.text(R_vals[0], W_vals[0], 'A', fontsize=12, ha='left', va='bottom')
    plt.text(R_vals[1], W_vals[1], 'B', fontsize=12, ha='left', va='bottom')
    plt.text(R_vals[2], W_vals[2], 'C', fontsize=12, ha='left', va='bottom')
    plt.text(R_vals[3], W_vals[3], 'D', fontsize=12, ha='left', va='bottom')
    plt.axhline( y = Max_Payload, color = 'r', linestyle = "--", linewidth = 2, label = "Max Payload")

    plt.xlabel('Range [nmi]')
    plt.xlim(left = -5)
    plt.ylabel('Payload Weight [lbf]')
    plt.title('Payload-Range Chart')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# Settings
Total_Mission_Fuel = 659.59 # lbf
Max_Payload = 3000           # lbf
W_e = 4839.33        # Operational empty eight: includes aircraft empty weight and pilot
CD_o = 0.04048
R_ferry = 3.646e+6           # ft (600 nmi)
e = 1.295
eta = 0.72
AR = 5
c = 0.6  # lbm/hp/hr

PlotPayloadRange(Total_Mission_Fuel, Max_Payload, W_e, CD_o, R_ferry, e, eta, AR, c)
