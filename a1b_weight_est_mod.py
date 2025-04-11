import numpy as np
from scipy.optimize import root_scalar

# Weight Estimate Calculator Modified for Assignment A1B

def CruiseFraction(R, c, eta, lift_to_drag):


    """ Calculates the weight fraction of cruise stage based on needed range
        R: The range the aircraft will travel during segment in ft
        c: Specific fuel consumption of engine in ft^-1
        eta: Estimated efficiency of propeller
        lift_to_drag: Estimated using wetted aspect ratio
        
        returns: weight fraction of cruise segment"""
    
    weight_frac = np.exp( (-R*SfcUnitConverter(c)) / (eta * lift_to_drag))
    #print("cruise",weight_frac)

    return weight_frac

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
    #print("loiter",weight_frac)

    return weight_frac

def NewWeight(crew_weight, payload_weight, empty_weight, fuel_weight_frac):

    """Calculates a new guess value for the takeoff weight of the aircraft.
        crew_weight: weight of the crew in lbf
        payload_weight: weight of the payload in lbf
        empty_weight: empty weight of aircraft in lbs
        fuel_weight_frac: calculated total fuel weight fraction"""

    weight = (crew_weight + payload_weight + empty_weight) / (1 - fuel_weight_frac)

    return weight

def SfcUnitConverter(sfc):
    """Converts a specific fuel consumption from lbm/hp/hr to ft^-1
        bsfc: Brake Specific fuel consumption in lbf/hp/hr"""
    sfc_convert = sfc/(1980000)

    return sfc_convert

def TaxiFraction(P_shaft, W, c_SL = 0.30):
    """Calculates W_1/W_o or the fuel burn from running the engine for 15 min @ 5% max
        P: The required thrust power of the engine in horsepower
        W: Initial gross weight of the aircraft
        prop_eta: The propeller efficiency
        c_sl: Predicted specific fuel consumption at sea level in lbf/hp*hr"""
    
    frac = 1 - (c_SL * P_shaft * 0.05 * 15 / (W * 60))
    #print("taxi",frac)

    return frac

def TakeoffFraction(P_shaft, W, c_SL):
    """Calculates W_2/W_1 or the fuel burn from running the engine at max power for 1 min
        P: The required thrust power of the engine in horsepower
        W: Gross weight of aircraft after taxi
        prop_eta: The propeller efficiency
        c_sl: Predicted specific fuel consumption at sea level in lbf/hp*hr"""
    
    frac = 1 - (c_SL * P_shaft / (W * 60))
    #print("takeoff",frac)

    return frac

def ZeroLiftDrag(S_ref, W_o, cf):
    """Calculates the zero lift drag of the aircraft based on the takeoff weight and wing area (applies for clean config only)
        S_ref: The wing area of the aircraft in ft^2
        W_o: The GTOW of the aircraft in lbs
        cf: The predicted skin friction factor of the aircraft"""

    # Roskam constants
    c = 1.0447
    d = 0.5326

    S_wet = (10**c) * (W_o**d)
    f = S_wet*cf

    Cdo = f/S_ref

    return Cdo

def LiftCoefficient(Cdo, AR, e):
    """Calculates the lift coefficient when the aircraft is in its clean configuration"""
    k = (np.pi * e* AR)**-1
    C_L = np.sqrt(Cdo/k)
    L_D = 0.94 * C_L / (Cdo + k*(C_L**2))

    return L_D


def ExponentialFuelFrac(W_i, S, V, eta, rho, AR, e, Cdo, c_t, R, n=50):
    """Discretized cruise fuel burn calculator.
    W_i: Initial weight at start of cruise segment (lbf)
    S: Reference wing area (ft^2)
    V: Cruise velocity (ft/s)
    rho: Air density (slugs/ft^3)
    AR: Aspect ratio
    e: Oswald efficiency
    Cdo: Zero-lift drag coefficient
    c_t: Specific fuel consumption (1/ft)
    R: Total cruise range (ft)
    n: Number of discretized segments

    Returns:
        Final cruise segment weight fraction (W_final / W_initial)
    """

    delta_R = R / n
    W = W_i
    k = 1 / (np.pi * e * AR)

    for _ in range(n):
        CL = (2 * W) / (rho * V**2 * S)
        CD = Cdo + k * CL**2
        LD = CL / CD
        dW_frac = np.exp(-c_t * delta_R / (eta * LD))
        W *= dW_frac

    return W / W_i

def StandardDensity(h):
    """Returns density in slugs/ft^3"""
    rho_SL = 2.377e-3  # slugs/ft^3
    delta = (1 - 6.875e-6 * h) ** 5.2561  # pressure ratio
    theta = 1 - 6.875e-6 * h  # temperature ratio
    sigma = delta / theta  # density ratio
    rho = sigma * rho_SL
    return rho

def ClimbFuelFrac(W_i, S, h_i, h_f, P_shaft, eta, CDo, e, AR, c, n = 50):

    g = 32.17 # ft/s^2
    CLmax = 1.6
    k = 1 / (np.pi * e * AR)
    eta_total = 0.7025
    P_hp = eta_total * P_shaft # available power
    P = P_hp * 550 # lbf*ft/s

    climb_h = h_f - h_i # change in altitude in ft

    delta_h = climb_h/n # height change during increment
    W = W_i

    current_altitude = h_i
    x_climb = 0 
    for _ in range(n):

        rho = StandardDensity(current_altitude)
        V = ((4/3) * (W/S)**2 * k / (CDo * rho**2))**(1/4)

        V_stall = np.sqrt((2 * W) / (rho * S * CLmax))

        if V_stall * 1.05 > V:
            V = V_stall

        T = P/V

        A = P * np.sqrt(0.5*rho*S) / W**(1.5)

        def f(CL):
            return k * CL**2 - (A/2) * CL**1.5 - CDo
        
        # Solve using root_scalar
        sol = root_scalar(f, bracket=[1e-6, 10.0], method='brentq')
        if not sol.converged:
            raise RuntimeError("Failed to converge to a CL solution.")
        CL = sol.root

        if CL > CLmax:
            CL = CLmax

        CD = CDo + k*(CL**2)
        q = 0.5 * rho * (V**2)
        D = q * S * CD

        del_he = delta_h + ((V**2) / (2*g))
        dW_frac = np.exp(-(c* del_he * T)/(eta* W *(1-(D/T))))
        current_altitude += delta_h

        P_s = V * (T-D) / W
        x_climb += delta_h * V / P_s
        W*= dW_frac

    

    return W/W_i,  x_climb




def solve_takeoff_weight_3(W_crew, W_payload, S, W_e, mission_segments,
                            tol=1e-6, max_iter=10000):

    takeoff_weight_guess = 10000  # Initial guess
    iterations = 0

    while iterations < max_iter:
        current_weight = takeoff_weight_guess

        # Aircraft constants
        AR = 5
        e = 1.192
        cf = 0.0150
        rho = 0.0022409  # or change by altitude
        Cdo_loiter = ZeroLiftDrag(S, takeoff_weight_guess, cf)
        L_D_loiter = LiftCoefficient(Cdo_loiter, AR, e)

        segment_fuel_fracs = []

        for seg in mission_segments:
            seg_type = seg["type"]

            if seg_type == "taxi":
                frac_taxi = TaxiFraction(seg["P_shaft"], current_weight, seg.get("c_SL", 0.3))
                segment_fuel_fracs.append(frac_taxi)
                current_weight *= frac_taxi

            elif seg_type == "takeoff":
                frac_takeoff = TakeoffFraction(seg["P_shaft"], current_weight, seg.get("c_SL", 0.6))
                segment_fuel_fracs.append(frac_takeoff)
                current_weight *= TakeoffFraction(seg["P_shaft"], current_weight, seg.get("c_SL", 0.6))

            elif seg_type == "climb":
                c_t = SfcUnitConverter(seg["c"])
                rho = StandardDensity(seg["h_i"])
                Cdo = ZeroLiftDrag(S, takeoff_weight_guess, cf)
                frac_climb, x_climb = ClimbFuelFrac(
                    W_i=current_weight,
                    S=S,
                    h_i=seg["h_i"],
                    h_f=seg["h_f"],
                    P_shaft=seg["P_shaft"],
                    eta=seg["eta"],
                    CDo=Cdo,
                    e=e,
                    AR=AR,
                    c=c_t
                )
                
                segment_fuel_fracs.append(frac_climb)
                current_weight *= frac_climb


            elif seg_type == "descend" or seg_type == "landing":
                segment_fuel_fracs.append(seg["fraction"])
                current_weight *= seg["fraction"]

            elif seg_type == "cruise":
                c_t = SfcUnitConverter(seg["c"])
                R = seg["R"] - x_climb if x_climb > 0 else seg["R"]
                cruise_frac = ExponentialFuelFrac(current_weight, S, seg["V"], seg["eta"], rho, AR, e, Cdo, c_t, R)

                segment_fuel_fracs.append(cruise_frac)
                current_weight *= cruise_frac
                x_climb = 0 # reset climb distance

            elif seg_type == "loiter":
                loiter_frac = LoiterFraction(seg["E"], seg["V"], seg["c"], seg["eta"], L_D_loiter)
                segment_fuel_fracs.append(loiter_frac)
                current_weight *= loiter_frac

            else:
                raise ValueError(f"Unknown segment type: {seg_type}")
        

        reserves = 1.06

        fuel_weight_frac = (1 - (current_weight / takeoff_weight_guess)) * reserves

        new_takeoff_weight = NewWeight(W_crew, W_payload, W_e, fuel_weight_frac)

        error = abs(new_takeoff_weight - takeoff_weight_guess) / new_takeoff_weight

        if 2 * error < tol:
            fuel_weight = fuel_weight_frac * new_takeoff_weight
            return (new_takeoff_weight, W_e/new_takeoff_weight, fuel_weight_frac,
                    W_e, fuel_weight, iterations, segment_fuel_fracs)

        fuel_weight = fuel_weight_frac * new_takeoff_weight
        takeoff_weight_guess = new_takeoff_weight

        iterations += 1

    raise ValueError(f"Did not converge within {max_iter} iterations. Last error: {error:.6f}")

################################################### CODE IMPLEMENTATION BEGINS HERE #######################################
# Aircraft Parameters

W_crew = 190 # lbs
W_payload = 3000 # lbs
S = 387.2 # ft^2
W_empty = 4881.4 # lbs
P_shaft = 750 # shp

# Define mission segments MUST BE SEQUENTIAL    
mission_segments = [
    {"type": "taxi", "P_shaft": P_shaft, "c_SL": 0.3},
    {"type": "takeoff", "P_shaft": P_shaft, "c_SL": 0.6},
    {"type": "climb", "h_i": 2000, "h_f": 3000, "P_shaft": P_shaft, "eta": 0.77, "c": 0.6},
    {"type": "cruise", "R": 151903, "c": 0.6, "eta": 0.72, "V": 213.255},
    {"type": "descend", "fraction": 0.999},
    {"type": "cruise", "R": 1033065, "c": 0.6, "eta": 0.72, "V": 180.446},
    {"type": "climb", "h_i": 2000, "h_f": 3000, "P_shaft": P_shaft, "eta": 0.77, "c": 0.6},
    {"type": "cruise", "R": 151903, "c": 0.6, "eta": 0.72, "V": 213.255},
    {"type": "loiter", "E": 0.5, "V": 180, "c": 0.45, "eta": 0.82},
    {"type": "descend", "fraction": 0.999},
    {"type": "landing", "fraction": 0.998}
]

results = solve_takeoff_weight_3(W_crew, W_payload, S, W_empty, mission_segments)

takeoff_weight, empty_weight_frac, fuel_weight_frac,\
empty_weight, fuel_weight, iterations, segment_fuel_fracs = results

print(f"Takeoff Weight Estimate after {iterations} iterations:")
print(f"  Takeoff Weight (W_naught): {takeoff_weight:.2f} lbs")
print(f"  Empty Weight Fraction: {empty_weight_frac:.4f}")
print(f"  Fuel Weight Fraction: {fuel_weight_frac:.4f}")
print(f"  Empty Weight: {empty_weight:.2f} lbs")
print(f"  Fuel Weight: {fuel_weight:.2f} lbs")

print("\nSegment Fuel Fractions:")
for i, frac in enumerate(segment_fuel_fracs, 1):
    print(f"  Segment {i}: Fuel burned = {1 - frac:.4f}")
