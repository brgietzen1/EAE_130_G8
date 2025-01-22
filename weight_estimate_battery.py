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
    


def solve_takeoff_weight_HCFA(crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, custom_segments, e=1e-6, max_iter=1000):
    """
    Iteratively solve for the takeoff weight of the aircraft using hydrocarbon fuels.
    
    Parameters:
    - crew_weight: Total crew weight (including baggage)
    - payload_weight: Total payload weight
    - A, C: Constants from Raymer tables for empty weight fraction
    - cruise_segments: List of tuples [(R, c, eta, lift_to_drag), ...] for cruise segments
    - loiter_segments: List of tuples [(E, V, c, eta, lift_to_drag), ...] for loiter segments
    - custom_segments: List of weight fractions for non-cruise/loiter segments
    - e: Convergence tolerance
    - max_iter: Maximum number of iterations
    
    Returns:
    - takeoff_weight: Converged takeoff weight
    - cruise_fuel_weight_HCFA: Total weight of the fuel needed for the cruise segments (in lbs)
    - fuel_fractions: List of fuel fractions for loiter segments
    """
    # Initial guess for takeoff weight
    takeoff_weight_guess = crew_weight + payload_weight + 1000  # Arbitrary guess
    iterations = 0
   

    while iterations < max_iter:
        loiter_fuel_fractions = []
        cruise_fuel_weight_HCFA = 0 # Track the fuel weight used in cruise segments
        # Step 1: Calculate the empty weight fraction
        empty_weight_frac = RaymerMethod(takeoff_weight_guess, A, C)

        # Step 2: Calculate the total fuel weight fraction
        fuel_weight_frac = 1.0  # Start with a neutral multiplier
        
        # Cruise segments
        for R, c, eta, lift_to_drag in cruise_segments:
            cruise_fuel_weight_HCFA += CruiseFraction(R, c, eta, lift_to_drag) * takeoff_weight_guess  # Add fuel used in this segment
            fuel_weight_frac *= CruiseFraction(R, c, eta, lift_to_drag)
        
        # Loiter segments
        for E, V, c, eta, lift_to_drag in loiter_segments:
            fuel_weight_frac *= LoiterFraction(E, V, c, eta, lift_to_drag)
            loiter_fuel_fractions.append(fuel_weight_frac)

        # Custom segments (These values come directly from the textbook)
        for frac in custom_segments:
            fuel_weight_frac *= frac

        # Step 3: Solve for the new takeoff weight
        new_takeoff_weight = NewWeight(crew_weight, payload_weight, fuel_weight_frac, empty_weight_frac)

        # Step 4: Calculate error
        error = abs(new_takeoff_weight - takeoff_weight_guess) / new_takeoff_weight

        # Step 5: Check convergence
        if error < e:
            return new_takeoff_weight, cruise_fuel_weight_HCFA, loiter_fuel_fractions

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



def EnergyFractions(custom_segments, loiter_fractions, hcfa_togw, cruise_fuel_weight_HCFA):
    """
    Calculates the energy fractions needed for estimating the weight of an electric aircraft.
    
    Parameters:
    - custom_segments: List of weight fractions for HCFA aircraft.
    - loiter_fractions: List of loiter weight fractions from the final HCFA solver iteration.
    - hcfa_togw: The final calculated takeoff weight of the HCFA aircraft.
    - fuel_weight_HCAS: The fuel weight needed for cruise (not the fuel weight fraction) in lbs for the HCFA craft. 
    
    Returns:
    - energy_frac: List of energy fractions for each flight segment.
    """
    energy_fractions = []
    
    # For custom segments, calculate energy fractions
    for frac in custom_segments:
        energy_fraction = (1 - frac) * (hcfa_togw / cruise_fuel_weight_HCFA)
        energy_fractions.append(energy_fraction)

    # For loiter segments, calculate energy fractions
    for loiter_frac in loiter_fractions:
        energy_fraction = (1 - loiter_frac) * (hcfa_togw / cruise_fuel_weight_HCFA)
        energy_fractions.append(energy_fraction)

    return energy_fractions
    

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

def solve_takeoff_weight_electric(crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, custom_segments, eta_battery, specific_energy, e=1e-6, max_iter=1000):
    """
    Iteratively solve for the takeoff weight of the electric aircraft using HCFA solver for energy fractions.
    
    Parameters:
    - crew_weight: Total crew weight (including baggage) in lbf
    - payload_weight: Total payload weight in lbf
    - A, C: Constants from Raymer tables for empty weight fraction
    - cruise_segments: List of tuples [(R, c, eta, lift_to_drag), ...] for cruise segments
    - loiter_segments: List of tuples [(E, V, c, eta, lift_to_drag), ...] for loiter segments
    - custom_segments: List of weight fractions for non-cruise/loiter segments
    - eta_battery: Battery efficiency (fractional)
    - specific_energy: Specific energy of the battery in Wh/kg
    - e: Convergence tolerance
    - max_iter: Maximum number of iterations
    
    Returns:
    - takeoff_weight: Converged takeoff weight for the electric aircraft in lbf
    - battery_mass: Mass of the battery required for the cruise segment in slugs
    - energy_fractions: List of energy fractions for each segment
    """
    
    # Initial guess for electric aircraft weight
    weight_guess = crew_weight + payload_weight + 1000  # Arbitrary initial guess
    iterations = 0
    battery_mass = 0  # Start with no battery mass

    while iterations < max_iter:

        energy_fractions = [] # create energy fraction list
        battery_mass = 0 # create variable for battery mass

        # Step 1: Solve for HCFA takeoff weight and cruise fuel weight
        hcfa_togw, cruise_fuel_weight_HCFA, loiter_fuel_fractions = solve_takeoff_weight_HCFA(
            crew_weight, payload_weight, A, C, cruise_segments, loiter_segments, custom_segments, e, max_iter
        )
        
        # Step 2: Calculate the empty weight fraction using the Raymer method
        W_e = RaymerMethod(weight_guess, A, C)

        # Step 3: Calculate the total energy fractions using the HCFA results
        energy_fractions = EnergyFractions(custom_segments, loiter_fuel_fractions, hcfa_togw, cruise_fuel_weight_HCFA)

        # Step 4: Calculate battery mass for cruise using the current weight guess
        battery_mass = 0
        for R, _, _, l_to_d in cruise_segments:
            battery_mass += CruiseBatteryMass(R, weight_guess, eta_battery, specific_energy, l_to_d)

        # Step 5: Solve for the new electric aircraft weight using the calculated values
        new_weight = NewElectricWeight(weight_guess, crew_weight, payload_weight, W_e, battery_mass, energy_fractions)

        # Step 6: Calculate error for convergence check
        error = abs(new_weight - weight_guess) / new_weight
        
        # Step 7: Check convergence
        if error < e:
            return new_weight, battery_mass, energy_fractions

        # Update guess and iterate
        weight_guess = new_weight
        iterations += 1

    # If the loop completes without convergence
    raise ValueError(f"Did not converge within {max_iter} iterations. Last error: {error:.6f}")

    
# Sample data for the electric aircraft solver

# Constants from Raymer tables for empty weight fraction calculation
A = 0.74   # Example constant from Raymer (this should be based on your aircraft design)
C = -0.03     # Example exponent from Raymer (this should be based on your aircraft design)

# Crew and payload weights in lbf (pounds force)
crew_weight = 240  # Weight of the crew and baggage in lbf
payload_weight = 3500  # Weight of the payload in lbf

# Cruise segments for the HCFA and electric aircraft
# Each tuple contains (R, c, eta, lift_to_drag)
cruise_segments = [
    (10000, 0.5, 0.85, 15),  # Example cruise segment: (Range in ft, Specific Fuel Consumption, Efficiency, Lift-to-Drag ratio)
    (15000, 0.6, 0.87, 16)   # Another example cruise segment
]

# Loiter segments for the HCFA and electric aircraft
# Each tuple contains (E, V, c, eta, lift_to_drag)
loiter_segments = [
    (1200, 150, 0.4, 0.9, 12),  # Example loiter segment: (Endurance in seconds, Velocity in ft/s, Specific Fuel Consumption, Efficiency, Lift-to-Drag ratio)
    (1500, 140, 0.45, 0.88, 14) # Another example loiter segment
]

# Custom segments (weight fractions) for other non-cruise/loiter segments
custom_segments = [.996,.995,.996,.998,.999,.998]  # Example custom weight fractions

# Battery efficiency for electric aircraft (fractional)
eta_battery = 0.95  # Example battery efficiency (95%)

# Specific energy of the battery in Wh/kg
specific_energy = 220  # Example specific energy (220 Wh/kg)

# Lift-to-drag ratio for the electric aircraft
lift_to_drag = 10  # Example lift-to-drag ratio (dimensionless)

# Convergence tolerance and maximum iterations for iterative solver
e = 1e-6
max_iter = 1000

# Solve for the takeoff weight, battery mass, and energy fractions for the electric aircraft
takeoff_weight, battery_mass, energy_fractions = solve_takeoff_weight_electric(
    crew_weight,
    payload_weight,
    A,
    C,
    cruise_segments,
    loiter_segments,
    custom_segments,
    eta_battery,
    specific_energy,
    e=e,
    max_iter=max_iter
)

# Print the results
print(f"Converged takeoff weight of the electric aircraft: {takeoff_weight:.2f} lbf")
print(f"Required battery mass for cruise segments: {battery_mass:.2f} slugs")
print("Energy fractions for each segment:")
for i, energy_frac in enumerate(energy_fractions):
    print(f"  Segment {i+1}: {energy_frac:.4f}")





