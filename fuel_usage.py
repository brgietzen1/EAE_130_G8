import numpy as np

def TaxiFraction(P, W, prop_eta = 0.77, c_SL = 0.60):
    """Calculates W_1/W_o or the fuel burn from running the engine for 15 min @ 5% max
        P: The required thrust power of the engine in horsepower
        W: Initial gross weight of the aircraft
        prop_eta: The propeller efficiency
        c_sl: Predicted specific fuel consumption at sea level in lbf/hp*hr"""
    
    P_shaft = P/prop_eta
    frac = 1 - (c_SL * P_shaft * 0.05 * 15 / (W * 60))

    return frac

def TakeoffFraction(P, W, prop_eta = 0.77, c_SL = 0.60):
    """Calculates W_2/W_1 or the fuel burn from running the engine at max power for 1 min
        P: The required thrust power of the engine in horsepower
        W: Gross weight of aircraft after taxi
        prop_eta: The propeller efficiency
        c_sl: Predicted specific fuel consumption at sea level in lbf/hp*hr"""
    
    P_shaft = P/prop_eta
    frac = 1 - (c_SL * P_shaft / (W * 60))

    return frac

def CruiseFraction(R, c, eta, lift_to_drag):


    """ Calculates the weight fraction of cruise stage based on needed range
        R: The range the aircraft will travel during segment in ft
        c: Specific fuel consumption of engine in ft^-1
        eta: Estimated efficiency of propeller
        lift_to_drag: Estimated using wetted aspect ratio
        
        returns: weight fraction of cruise segment"""
    
    weight_frac = np.exp( (-R*SfcUnitConverter(c)) / (eta * lift_to_drag))

    return weight_frac

def SfcUnitConverter(sfc):
    """Converts a specific fuel consumption from lbm/hp/hr to ft^-1
        bsfc: Brake Specific fuel consumption in lbf/hp/hr"""
    sfc_convert = sfc/(1980000)

    return sfc_convert

def LoiterFraction(E, V, c, eta, lift_to_drag):

    """ Calculates the weight fraction of the loiter stage based on needed endurance
        E: Time aircraft will loiter during segment in seconds
        V: Expected velocity of craft during loiter period in ft/s
        c: Specific fuel consumption of engine in ft^-1
        eta: Estimated efficiency of propeller
        lift_to_drag: Estimated using wetted aspect ratio
        
        returns: weight fraction of loiter segment"""
    
    E_seconds = E*3600

    weight_frac = np.exp( (-E_seconds*V*SfcUnitConverter(c)) / (eta*lift_to_drag))

    return weight_frac

W_TO = 8864
P_TO = 620
lift_to_drag = 11.648

taxi_frac = TaxiFraction(P_TO, W_TO, prop_eta = 0.77, c_SL = 0.60)
fuel_weight_taxi = W_TO - (taxi_frac*W_TO)
print("Fuel used in taxi:",fuel_weight_taxi)

W1 = (taxi_frac*W_TO)
takeoff_frac = TakeoffFraction(P_TO, W1, prop_eta = 0.77, c_SL = 0.60)
fuel_weight_takeoff = W1 - takeoff_frac*W1

print("Fuel used in takeoff:",fuel_weight_takeoff)

W2 = takeoff_frac*W1
cruise_frac1 = CruiseFraction(151903, 0.6, 0.82, lift_to_drag)
fuel_weight_cruise1 = W2 - cruise_frac1*W2
print("Fuel used in cruise 1:",fuel_weight_cruise1)

W3 = cruise_frac1*W2
cruise_frac2 = CruiseFraction(1033065, 0.6, 0.82, lift_to_drag)
fuel_weight_cruise2 = W3 - cruise_frac2*W3
print("Fuel used in cruise 2:",fuel_weight_cruise2)

W4 = cruise_frac2*W3
cruise_frac3 = CruiseFraction(151903, 0.6, 0.82, lift_to_drag)
fuel_weight_cruise3 = W4 - cruise_frac3*W4
print("Fuel used in cruise 3:",fuel_weight_cruise3)


W5 = cruise_frac3*W4
loiter_frac = LoiterFraction(0.5, 180, 0.6, 0.72, lift_to_drag)
fuel_weight_loiter = W5- loiter_frac*W5
print("Fuel used in loiter:",fuel_weight_loiter)