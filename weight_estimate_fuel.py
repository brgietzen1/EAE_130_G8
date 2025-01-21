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
def solve_takeoff_weight(crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, custom_segments, e = 1e-6, max_iter=1000):
    """
    Iteratively solve for the takeoff weight of the aircraft.
    
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

# Example Usage
crew_weight = 240  # lbf
payload_weight = 3500  # lbf
A = 0.74              # Raymer constant for "Agricultural aircraft" (Table 3.1)
C = -0.03             # Raymer constant for "Agricultural aircraft" (Table 3.1)
tol = 1e-6            # Convergence tolerance

# Define flight segments
cruise_segments = [
    (600000, 0.25, 0.8, 10),  # Range in ft, specific fuel consumption in lbm/hp/hr, propeller efficiency, L/D
]
loiter_segments = [
    (0.5, 300, 0.40, 0.8, 10),  # Endurance in hrs, velocity in ft/s, specific fuel consumption lbm/hp/hr, propeller efficiency, L/D
]
other_segments = [0.97, 0.995, 0.98]  # Estimated weight fractions for non-cruise/loiter phases

# Solve for takeoff weight
takeoff_weight, iterations = solve_takeoff_weight(
    crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, other_segments, tol
)

print(f"Takeoff Weight: {takeoff_weight:.2f} lbf after {iterations} iterations.")


