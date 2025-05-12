import pandas as pd
import matplotlib.pyplot as plt

# Load CSVs for each flap configuration
df_clean = pd.read_csv(r"C:\Users\brgie\Documents\EAE 130\AVL\avl3.40b\Avl\runs\run_cases\polar_0f.csv")
df_takeoff = pd.read_csv(r"C:\Users\brgie\Documents\EAE 130\AVL\avl3.40b\Avl\runs\run_cases\polar_10f.csv")
df_landing = pd.read_csv(r"C:\Users\brgie\Documents\EAE 130\AVL\avl3.40b\Avl\runs\run_cases\polar_30f.csv")

# Group into a dictionary with colors
configs = {
    'Clean (0° Flap)': {'df': df_clean, 'color': 'blue'},
    'Takeoff (10° Flap)': {'df': df_takeoff, 'color': 'green'},
    'Landing (30° Flap)': {'df': df_landing, 'color': 'red'}
}

# Plot 1: Drag Polar (Cd_total vs CL)
plt.figure()
for label, cfg in configs.items():
    df = cfg['df']
    plt.plot(df['Cdtotal'], df['CL'], label=label, color=cfg['color'])
plt.xlabel('$C_D$')
plt.ylabel('$C_L$')
plt.title('Ceres-100 Drag Polar (Empty)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot 2: Induced Drag vs CL
"""plt.figure()
for label, cfg in configs.items():
    df = cfg['df']
    plt.plot(df['CL'], df['Cdind'], label=label, color=cfg['color'])
plt.xlabel('CL (Lift Coefficient)')
plt.ylabel('$C_D_}$')
plt.title('Induced Drag vs CL')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()"""

# Plot 3: Oswald Efficiency vs CL
plt.figure()
for label, cfg in configs.items():
    df = cfg['df']
    plt.plot(df['CL'], df['Oswald_e'], label=label, color=cfg['color'])
plt.xlabel('$C_L$')
plt.ylabel('Span Efficiency (e)')
plt.title('Span Efficiency vs $C_L$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
