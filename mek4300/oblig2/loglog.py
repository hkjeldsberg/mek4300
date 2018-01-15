import numpy as np
import matplotlib.pyplot as plt

# Parametres
rho = 1000.                # kg / m^3
U = 1.                     # m/s
mu = 8.9*10**(-4)          # Pa s
Re = np.linspace(1000,10000) # Dimensionless
L = mu*Re / (rho*U)        # m         

D = np.sqrt(rho*mu*np.pi*L*U**3* ( 2/np.pi - 0.5))
C_D = D / (0.5*rho*U**2*L)

# Create log-log plot
logCD = np.log(C_D)
logRe = np.log(Re)
plt.plot(Re, logCD)
plt.xlabel("log ( Reynolds number )")
plt.ylabel("log ( Drag coefficient )")
plt.show()

