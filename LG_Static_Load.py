import numpy as np

# Distance from nose landing gear to AFT cg
Na = 17.651
# Distance between nose and rear landing gear
B = 19.275
# Distance between AFT cg and rear landing gear 
Ma = 2.259
# Distance between FWD cg and rear landing gear 
Mf = 3.819
# Height of the landing gear to cg 
H = 6.838

W = 8893.5
W_landing = 7114.8
FWD = 17.382
AFT = 18.575
g = 32.17 #ft/s^2
V = 168.781 #ft/s

Max_Main_Static_Load = W * (Na / B)
Max_Nose_Static_Load = W * (Mf / B)
Min_Nose_Static_Load = W * (Ma / B)
Dynamic_Breaking_Load = 0.31 * (H / B) * W
KE = 0.5 * (W_landing / g) * V**2

print(f"Max Main Static Load is {Max_Main_Static_Load:.2f}")
print(f"Max Nose Static Load is {Max_Nose_Static_Load:.2f}")
print(f"Min Nose Static Load is {Min_Nose_Static_Load:.2f}")
print(f"Dynamic Breaking Load is {Dynamic_Breaking_Load:.2f}")
print(f"Kinetic Energy during landing is {KE:.2f} J")





