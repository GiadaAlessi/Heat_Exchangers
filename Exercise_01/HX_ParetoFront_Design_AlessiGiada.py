import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Initial Data
T1i = 95 + 273.15
T2i = 15 + 273.15
v1i = 1
v2i = 5
L = 0.5
a = 0.003
b = 0.05
c = 0.01
ew = 0.001
kw = 43

# Outlet Temp. Initial  Guess
T1og = 60 + 273.15
T2og = 40 + 273.15

# Fins parameters
kf = 43
heights = np.linspace(0.001, 0.01, 100)
thicknesses = np.linspace(0.0001, 0.004, 100)
num_fins = np.arange(10, 100, 90)
Lf = L

# Outlet Temp. Initial  Guess
T1og = 60 + 273.15
T2og = 40 + 273.15

# Iteration Process
tol = 1e-3
max_iter = 100

# Number of Configurations and Result's Array
n_configurations = 500
results = []

# Pareto Front 
for _ in range(n_configurations):

    # Casual Selection of Parameters
    hf = np.random.choice(heights)
    ef = np.random.choice(thicknesses)
    Nf = np.random.choice(num_fins)
    
    # Select only valid configurations
    if Nf * ef >= b:
        continue

    # Geometry Parameters
    S1 = a * b
    S2fix = c * b
    S2 = c * b - (Nf * ef * hf)
    Af = 2 * Nf * hf * L + Nf * ef * L
    Aoeff = b * L - (Nf * ef * L)
    Ao = Aoeff + Af
    Ai = b * L
    Dh1 = 4 * S1 / (2 * a + 2 * b)
    Dh2 = 4 * S2 / (2 * c + 2 * b + Nf * hf * 2)

    for iteration in range(max_iter):

        # Average Temperatures
        T1m = (T1i + T1og) / 2
        T2m = (T2i + T2og) / 2

        # Physical Properties
        rho1 = PropsSI('D', 'T', T1m, 'P', 101325, 'Water')
        cp1 = PropsSI('C', 'T', T1m, 'P', 101325, 'Water')
        mu1 = PropsSI('V', 'T', T1m, 'P', 101325, 'Water')
        muw = PropsSI('V', 'T', T1i, 'P', 101325, 'Water')
        k1 = PropsSI('L', 'T', T1m, 'P', 101325, 'Water')
        rho2 = PropsSI('D', 'T', T2m, 'P', 101325, 'Air')
        cp2 = PropsSI('C', 'T', T2m, 'P', 101325, 'Air')
        mu2 = PropsSI('V', 'T', T2m, 'P', 101325, 'Air')
        k2 = PropsSI('L', 'T', T2m, 'P', 101325, 'Air')

        # Mass Flow Rate
        mf1 = rho1 * v1i * S1
        mf2 = rho2 * v2i * S2fix

        # Capacities
        C1 = mf1 * cp1
        C2 = mf2 * cp2
        Cmin = min(C1, C2)
        Cmax = max(C1, C2)
        Z = Cmin / Cmax

        # Heat Transfer Coefficients
        Re1 = rho1 * v1i * Dh1 / mu1
        Re2 = rho2 * v2i * Dh2 / mu2
        Pr1 = mu1 * cp1 / k1
        Pr2 = mu2 * cp2 / k2
        v1o = mf1/(rho1*S1)
        v2o = mf2/(rho2*S2)
        vm1 = (v1i+v1o)/2
        vm2 = (v2i+v2o)/2
        m1 = 0.8
        m2 = 0.8
        n1 = 0.33
        n2 = 0.4
        c1 = 0.027
        c2 = 0.023
        K1 = (mu1/muw)**0.14
        K2 = 1
        Nu1 = c1 * Re1**m1 * Pr1**n1 * K1
        Nu2 = c2 * Re2**m2 * Pr2**n2 * K2
        a2 = Nu2 * k2 / Dh2
        a1 = Nu1 * k1 / Dh1

        # Skip invalid configurations
        if kf <= 0 or ef <= 0 or a2 <= 0:
            continue


        # Fin Efficiencies
        m = np.sqrt((2*a2)/(kf*ef))
        etaf = np.tanh(m * hf)/(m * hf)
        etao = 1- (Af / Ao)*(1 - etaf)

        # Overall Heat Transfer Coefficient
        Uo = (1/a1+ew/kw+1/(a2*etao))**-1

        # Epsilon - NTU
        NTU = Uo * Ao / Cmin
        epsilon = (1-np.exp(-NTU*(1-Z)))/(1-Z*np.exp(-NTU*(1-Z)))
        
        # Heat Exchanged
        Q = epsilon * Cmin * (T1i - T2i)

        # Evaluate Outlet Temperatures
        T1o = T1i - Q / C1
        T2o = T2i + Q / C2

        # Pressure Drops
        f1 = 0.0625 / ((np.log10(5.74/(Re1**0.9)))**2)
        f2 = 0.0625 / ((np.log10(5.74/(Re2**0.9)))**2)
        tau1 = f1 * rho1 * vm1**2 / 2
        tau2 = f2 * rho2 * vm2**2 / 2
        Dp1 = (mf1 * (v1o-v1i) + tau1 * Ai) / S1
        Dp2 = (mf2 * (v2o-v2i) + tau2 * Ao) / S2

        # Convergence check
        if abs(T1o - T1og) < tol and abs(T2o - T2og) < tol:
            break

        # Update Guess Outlet Temperature, if necessary
        T1og = T1o
        T2og = T2o
    
    # Filter Results
    if (Dp2) <= 200:
        results.append([Q, Dp2])

# Results
results = np.array(results)

# Pareto Front Graph
plt.figure(figsize=(10, 6))
plt.scatter(results[:, 0], results[:, 1], c='blue', alpha=0.6, label="Simulated Configurations")
plt.title("Pareto Front: Heat Exchanged vs Pressure Drop")
plt.xlabel("Q [W]")
plt.ylabel("ΔP air side [Pa]")
plt.grid(True)
plt.legend()
plt.show()


# Pareto front identification
def find_pareto_front(results):
    pareto_front = []
    for i, point in enumerate(results):
        dominated = False
        for other_point in results:
            # Check if the current point is dominated by another point
            if other_point[0] >= point[0] and other_point[1] <= point[1] and not np.array_equal(point, other_point):
                dominated = True
                break
        if not dominated:
            pareto_front.append(point)
    return np.array(pareto_front)

# Identify Pareto front
pareto_front = find_pareto_front(results)

# Sort Pareto front by heat exchanged
pareto_front = pareto_front[pareto_front[:, 0].argsort()]

# Pareto Front Graph
plt.figure(figsize=(10, 6))
plt.scatter(results[:, 0], results[:, 1], c='blue', alpha=0.6, label="Simulated Configurations")
plt.scatter(pareto_front[:, 0], pareto_front[:, 1], c='red', label="Pareto Front")
plt.title("Pareto Front: Heat Exchanged vs Pressure Drop")
plt.xlabel("Q [W]")
plt.ylabel("Δp Air side [Pa]")
plt.grid(True)
plt.legend()
plt.show()
