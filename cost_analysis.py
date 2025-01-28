# EAE130 Cost Analysis

import numpy as np

# (A) RDT&E COSTS
# Parameters Variables
#empty weight (lb) - !from weight calc
W_e = empty_weight_frac * new_takeoff_weight
#production quantity
Q =
#maximum velocity (knots)
V =
#engine maximum thrust (lb)
T_max =
#max mach number
M_max =
#turbine inlet temperature (Rankine)
T_turbine_inlet =
#no. of total engines per aircraft
n_engines = 
#total production quantity times number of engines per aircraft
N_eng_total = Q * n_engines
# Number of flight test aircraft
FTA =
#Engineering hourly rate (wrap rate - salary plus all other costs like benefits and administrative costs - typically the salary is half the wrap)
R_engineer = 2.576 * 2025 - 5058
#Tooling hourly rate
R_tooling = 2.883 * 2025 - 5666
#Quality control hourly rate 
R_qc = 2.60 * 2025 - 5112
#Manufacturing hourly rate 
R_manufacturing = 2.316 * 2025 - 4552
#Airline Factor
AF =
#Route Factor
K = 
#Max Takeoff Weight (lb) - !need from weight calc
MTOW = 
#Block time in hrs (Total time aircraft is in use for mission - from wheel block removal to wheel block placement)
tb = 
#current year (year our technology level is considered at - our date of entry into service is 2035 so we should use 2035)
t_year = 2025
#base year (ideally use 1989 like the text book)
b_year = 1989
#Cost escalation factor (t_CEF/b_CEF)
CEF = (5.17053 + 0.104981 * (t_year - 2006))/(5.17053 + 0.104981 * (b_year - 2006))
#Fuel Weight (lb) - !from weight calc
W_f = fuel_weight_frac * new_takeoff_weight
#Price per gallon of fuel 
P_f = 
#fuel density (lbs/gal)
rho_f = 
#Battery weight (lb) - !from weight calc
W_b = battery_mass
#Price of electricity in $/kWh
P_elec = 
#Specific energy of battery in Wh/kg (0 is no batteries)
e_elec = 
#Range in nmi
R = 
#Airframe weight (Empty weight minus engine weight, battery weight, and motor weight) - !from weight calc
W_A = W_e - 
#Maintenance labor weigth in USD/hr for the year of interest
RL = 
#Maximum engine thrust in lbs
T0 = 
#Total takeoff shaft horsepower (all engines added up)
SHP_TO = 
# no. of turboprop engines
n_engines_turbo = 
# no. of hours between engine overhauls (usually between 3000 and 5000)
H_em = 
#Aircraft residual value factor (estimated over lifetime)
K_depreciation = 
#Number of years the aircraft is used
n = 
#hull insurance rate, usually assumed to be 2%
IR_a = 0.02 

# Define the functions for various cost calculations based on the provided formulas

#Engineering Hours
def engineering_hours(W_e, Q):
   engineering_hours = 4.86 * (W_e ** 0.777) * (V ** 0.894) * (Q ** 0.163)
   return engineering_hours


#Toolings Hours
def tooling_hours(W_e, V, Q):
   tooling_hours = 5.99 * (W_e ** 0.777) * (V ** 0.696) * (Q ** 0.263)
   return tooling_hours


#Manufacturing Hours
def manufacturing_hours(W_e, V, Q):
   manufacturing_hours = 7.37 * (W_e ** 0.82) * (V ** 0.484) * (Q ** 0.641)
   return manufacturing_hours


#QC Hours
def qc_hours(manufacturing_hours):
   qc_hours = 0.133 * manufacturing_hours
   return qc_hours

#Total RDT&E Cost
def RDTE_cost(engineering_hours, R_engineer, tooling_hours, R_tooling, manufacturing_hours, R_manufacturing, qc_hours, R_qc):
    RDTE_cost = engineering_hours * R_engineer + tooling_hours * R_tooling + manufacturing_hours * R_manufacturing + qc_hours * R_qc
    return RDTE_cost

# (B) Flyaway/Production

#Development Support Cost
def development_support_cost(W_e, V):
   development_support_cost = 45.42 * (W_e ** 0.630) * (V ** 1.3)
   return development_support_cost


#Flight Test Cost
def flight_test_cost(W_e, V, FTA):
   flight_test_cost = 1243.03 * (W_e ** 0.325) * (V ** 0.822) * (FTA ** 1.21)
   return flight_test_cost


#Manufacturing Materials Cost
def manufacturing_materials_cost(W_e, V, Q):
   manufacturing_materials_cost = 11.0 * (W_e ** 0.921) * (V ** 0.621) * (Q ** 0.799)
   return manufacturing_materials_cost


#Engineering Production Cost per engine
def engine_production_cost(T_max, M_max, T_turbine_inlet):
   engine_production_cost = 1548 * (0.043 * T_max + 243.25 * M_max) + 0.969 * T_turbine_inlet - 2228
   return engine_production_cost

#Total flyaway cost (neglects avionics for now)
def C_flyaway(development_support_cost, flight_test_cost, manufacturing_materials_cost, engine_production_cost, N_eng_total):
    C_flyaway = development_support_cost + flight_test_cost + manufacturing_materials_cost + engine_production_cost * N_eng_total
    return C_flyaway

#Total Aircraft Cost (to recoup RDT&E and Production)
def C_aircraft(RDTE_cost, C_flyaway):
    C_aircraft = RDTE_cost + C_flyaway
    return C_aircraft

#Engine cost per plane
def C_engine(SHP_TO, CEF):
    C_engine = 10^(2.5262 + 0.9465 * np.log10(SHP_TO)) * CEF
    return C_engine

#Airframe cost per plane (used for later cost calculations)
def C_airframe(C_aircraft, C_engine):
    C_airframe = C_aircraft - C_engine
    return C_airframe

# (C) DOC: Direct Operating Costs
# DOC = COC + FOC

# Crew
def C_crew(AF, K, MTOW, tb, CEF):
    C_crew = AF * (K * (MTOW)^0.4 * tb) * CEF
    return C_crew

# Attendants (Not applicable)

# Fuel
def C_fuel(W_f, P_f, rho_f):
    C_fuel = 1.02 * W_f * (P_f/rho_f)
    return C_fuel

# For hybrid electric propulsion
def C_electric(W_b, P_elec, e_elec):
    C_electric = 1.05 * W_b * P_elec * e_elec
    return C_electric

# Oil (neglected for now)
'''def C_oil(W_f, tb, P_oil, rho_oil)
W_oil = 0.0125 * W_f * (tb/100)
C_oil = 1.02 * W_oil * (P_oil/rho_oil)'''

# Landing Fees
def C_airport(MTOW, CEF):
    C_airport = 1.5 *(MTOW/1000) * CEF
    return C_airport

# Navigation Fees
def C_navigation(CEF, R, tb, MTOW):
    C_navigation = 0.5 * CEF * ((1.852 * R)/(tb)) * ((0.00045359237 * MTOW)/(50))^0.5
    return C_navigation

# Airframe Maintenance
def C_AirMain(W_A, RL, CEF, C_airframe, tb):
    C_ML = 1.03 * (3+ (0.067 * W_A)/(1000))* RL
    C_MM = 1.03 * (30 * CEF) + 0.79 * 10^(-5) * C_airframe
    C_AirMain = (C_ML + C_MM) * tb
    return C_AirMain

# Engine Maintenance labor cost
def C_ML_engine(T0, tb, RL):
    C_ML_engine = (0.645 + (0.05 * T0/10^4))* (0.566 + 0.434/tb)* RL
    return C_ML_engine

#Engine maintenance material cost
def C_MM_engine(T0, tb, CEF):
    C_MM_engine = ((25 + (18*T0/10^4)) * (0.62 + 0.38/tb)) * (CEF)
    return C_MM_engine

# for turboprop engines
def C_ML_turbo(SHP_TO, n_engines_turbo, H_em, RL):
    C_ML_turbo = 1.03 * 1.3 *(0.4956 + 0.0532 * (SHP_TO/n_engines_turbo)/(1000)* ((1100)/(H_em))+ 0.1) * RL
    return C_ML_turbo

#Total engine maintenance
def C_EngMain(n_engines, n_engines_turbo, C_ML, C_MM, C_ML_turbo, tb):
    C_EngMain = ((n_engines - n_engines_turbo) * C_ML + n_engines * C_MM + n_engines_turbo * C_ML_turbo) * tb
    return C_EngMain

#Electric motor cost
def C_motors(motor_hp):
    C_motors = 150*motor_hp
    return C_motors

#Battery cost
def C_battery(battery_kWh):
    C_battery = 520 * battery_kWh
    return C_battery

#Annual Utilization
def U_annual(tb):
    U_annual = 1.5 * 10^3 * (3.4546 * tb + 2.994 - (12.289 * tb^2 - 5.6626 * tb + 8.964)^0.5)
    return U_annual

# Insurance
def C_insurance(tb, IR_a, C_aircraft, U_annual):
    C_insurance = ((IR_a * C_aircraft)/(U_annual)) * tb
    return C_insurance

# Financing
def C_financing(DOC):
    C_financing = 0.07 * DOC
    return C_financing

# Depreciation
def C_depreciation(C_unit, K_depreciation, tb, n, U_annual):
    C_depreciation = (C_unit * (1- K_depreciation) * tb)/(n * U_annual)
    return C_depreciation

#Direct Operating Costs
def DOC(C_crew, C_fuel, C_electric, C_airport, C_navigation, C_AirMain, C_EngMain, C_battery, C_motor, C_insurance, C_depreciation):
    DOC = C_crew + C_fuel + C_electric + C_airport + C_navigation + C_AirMain + C_EngMain + C_battery + C_motor + C_insurance + C_depreciation
    return DOC

# Registration taxes
def C_registration(MTOW, DOC):
    C_registration = (0.001 + 10^(-8) * MTOW) * DOC
    return C_registration

def DOC_total(DOC, C_registration):
    DOC_total = DOC + C_registration
    return DOC_total