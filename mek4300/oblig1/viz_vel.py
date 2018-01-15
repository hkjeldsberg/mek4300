import numpy as np
from matplotlib.pyplot import plot, show, xlabel, title, ylabel

# Parameters
g = 9.81
h1 = h2 = 0.01
nu1 = 10.**(-6)
nu2 = 10.*nu1
rho2rho1 = 0.01
nu2nu1 = 10.
theta = np.pi/6

# Velocity profiles
def u_1(y):
    frac = (rho2rho1*nu2nu1) / (h2 + h1*rho2rho1*nu2nu1)
    return g*np.sin(theta)*( (h1**2 - y**2)/(2*nu1) + frac * ( h2**2/(2*nu2) - h1**2/(2*nu1))*(h1 + y))  

def u_2(y):
    frac1 = (rho2rho1*nu2nu1) / (h2 + h1*rho2rho1*nu2nu1)
    frac2 = 1. / (h2 + h1*rho2rho1*nu2nu1)
    return g*np.sin(theta)*( h1**2/(2*nu1) - y**2/(2*nu2) + (frac1*h1 + frac2*y)*( h2**2/(2*nu2) - h1**2/(2*nu1)))  

# Visualize
y1 = np.linspace(-h1, 0)
y2 = np.linspace(0, h2)

plot(y1, u_1(y1))
xlabel("Position [m]"); ylabel("Velocity [m/s]")
title("Velocity profile of lower fluid")
show()
plot(y2, u_2(y2))
xlabel("Position [m]"); ylabel("Velocity [m/s]")
title("Velocity profile of upper fluid")
show()
