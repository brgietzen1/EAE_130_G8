import numpy as np

import numpy as np

# Weight Estimate Calculator for Fueled Propeller Aircraft


def RaymerMethod(W_naught,A,C):
    
    """Uses existing regression analyses to calculate the empty weight of the aircraft
        values of needed constants can be found in the Raymer textbook 
        
        W_naught: Initial gross weight input in lbf
        A: Constant from Raymer tables
        C: Constant from Raymer tables
        
        returns: empty weight fraction of aircraft"""
    
    empty_weight_frac = A*(W_naught**C)

    return empty_weight_frac

def CruiseFraction(R, c, eta, lift_to_drag):


    """ Calculates the weight fraction of cruise stage based on needed range
        R: The range the aircraft will travel during segment in ft
        c: Specific fuel consumption of engine in ft^-1
        eta: Estimated efficiency of propeller
        lift_to_drag: Estimated using wetted aspect ratio
        
        returns: weight fraction of cruise segment"""
    
    weight_frac = np.exp( (-R*SfcUnitConverter(c)) / (eta * lift_to_drag))

    return weight_frac

def LoiterFraction(E, V, c, eta, lift_to_drag):

    """ Calculates the weight fraction of the loiter stage based on needed endurance
        E: Time aircraft will loiter during segment in seconds
        V: Expected velocity of craft during loiter period in ft/s
        c: Specific fuel consumption of engine in ft^-1
        eta: Estimated efficiency of propeller
        lift_to_drag: Estimated using wetted aspect ratio
        
        returns: weight fraction of loiter segment"""
    

    weight_frac = np.exp( (-E*V*SfcUnitConverter(c)) / (eta*lift_to_drag))

    return weight_frac

def NewWeight(crew_weight, payload_weight, fuel_weight_frac, empty_weight_frac):

    """Calculates a new guess value for the takeoff weight of the aircraft.
        crew_weight: weight of the crew in lbf
        payload_weight: weight of the payload in lbf
        fuel_weight_frac: calculated total fuel weight fraction
        empty_weight_fraction: calculated empty weight fraction"""

    weight = (crew_weight + payload_weight) / (1 - fuel_weight_frac - empty_weight_frac)

    return weight

def SfcUnitConverter(sfc):
    """Converts a specific fuel consumption from lbm/hp/hr to ft^-1
        sfc: Specific fuel consumption in lbm/hp/hr"""
    sfc_convert = sfc/(550*32.174)

    return sfc_convert
    


# Iterative Solver for Takeoff Weight
def solve_takeoff_weight_HCFA(crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, custom_segments, e = 1e-6, max_iter=1000):
    """
    Iteratively solve for the takeoff weight of the aircraft using hydrocarbon fuels. 
    You will need to use this so that you can compute the energy ratios when estimating the weight
    of an electric or hybrid electric aircraft. 
    
    Parameters:
    - crew_weight: Total crew weight (including baggage)
    - payload_weight: Total payload weight
    - A, C: Constants from Raymer tables for empty weight fraction
    - cruise_segments: List of tuples [(R, c, eta, lift_to_drag), ...] for cruise segments
    - loiter_segments: List of tuples [(E, V, c, eta, lift_to_drag), ...] for loiter segments
    - other_segments: List of weight fractions [guess_1, guess_2, ...] for non-cruise/loiter segments
    - tol: Convergence tolerance
    - max_iter: Maximum number of iterations
    
    Returns:
    - takeoff_weight: Converged takeoff weight
    - iterations: Number of iterations used
    """
    # Initial guess for takeoff weight
    takeoff_weight_guess = crew_weight + payload_weight + 1000  # Arbitrary guess
    iterations = 0

    while iterations < max_iter:
        # Step 1: Calculate the empty weight fraction
        empty_weight_frac = RaymerMethod(takeoff_weight_guess, A, C)

        # Step 2: Calculate the total fuel weight fraction
        fuel_weight_frac = 1.0  # Start with a neutral multiplier
        
        # Cruise segments
        for R, c, eta, lift_to_drag in cruise_segments:
            fuel_weight_frac *= CruiseFraction(R, c, eta, lift_to_drag)

            

        # Loiter segments
        for E, V, c, eta, lift_to_drag in loiter_segments:
            fuel_weight_frac *= LoiterFraction(E, V, c, eta, lift_to_drag)
        
        
        # Custom segments (These values come directly metabook/textbook)
        for frac in custom_segments:
            fuel_weight_frac *= frac

        # Step 3: Solve for the new takeoff weight
        new_takeoff_weight = NewWeight(crew_weight, payload_weight, fuel_weight_frac, empty_weight_frac)

        # Step 4: Calculate error
        error = abs(new_takeoff_weight - takeoff_weight_guess) / new_takeoff_weight

        # Step 5: Check convergence
        if error < e:
            return new_takeoff_weight, iterations 

        # Update guess and iterate
        takeoff_weight_guess = new_takeoff_weight
        iterations += 1

    # If the loop completes without convergence
    raise ValueError(f"Did not converge within {max_iter} iterations. Last error: {error:.6f}")

def SpecificEnergyConverter(e):

    """Converts an input specific energy in Wh/kg to lbf*ft/slug
        e: Specific energy of battery in Wh/kg
        
        returns: specific energy of battery in lbf*ft/slug"""
    
    e = e*38750

    return e
    
    
def CruiseBatteryMass(R, weight_guess, eta_battery, e, lift_to_drag):


    """ Calculates the necessary battery mass for an electric aircraft during cruise 
        R: Range of the aircraft during cruise in ft. 
        weight_guess: The togw weight estimate of the electric aircaft during the current iteration
        eta_batttery: Battery efficiency
        e: Specific energy of the battery in Wh/kg
        lift_to_drag: Estimated lift to drag ratio based on historical aircraft
        
        returns: Mass of battery needed for cruise segment in slugs"""
    
    m_battery = (R*weight_guess) / (eta_battery*SpecificEnergyConverter(e)*lift_to_drag)

    return m_battery



def EnergyFractions(custom_segments, loiter_weight_segments, hcfa_togw, fuel_weight_HCAS):

    """Calculates the energy fractions needed for estimating the weight of an electric aircraft. There should be one energy fraction
    for each segment of the flight excluding cruise segments. 
        custom_segments: The weight fuel fraction for HCFA aircraft which can be found in Table 2.2 of the metabook
        (This should be the same as the custom segments input into the HCFA solver)
        loiter_weight segments: A list of loiter weight segments computed by the last iteration of the HCAS solver
        that can be used as fuel fractions when no value can be foun in Table 2.2. 
        hcas_togw: The final calculated togw computed by the HCFA solver
        fuel_weight_HCAS: The fuel weight (not the fuel weight fraction) in lbs of the fuel needed for cruise
        
        
        returns: a list of the estimated energy fractions for each flight segment"""
    

    #energy_frac = (1 - flight_segment)(hcfa_togw/fuel_weight_HCAS)

    #return energy_frac
    

def NewElectricWeight(weight_guess, W_crew, W_payload, W_e, m_battery, energy_fractions):

    """Solves for the new weight of an electric aircraft using an initial guess weight
        weight_guess: The current guessed weight of the aircraft for current iteration
        W_crew: The weight of the crew in lbf
        W_payload: The weight of the payload in lbf
        W_e: The estimated empty mass fraction of the vehicle
        m_battery: The estimated mass of the battery in slugs
        energy_fractions: list of calculated energy fractions
        
        returns: the new weight of the electric aircraft in lbf"""
    
    new_weight = (W_crew*W_payload) / (1- W_e - ((m_battery*32.17/weight_guess)*(1 + np.sum(energy_fractions))))


    return new_weight
    




