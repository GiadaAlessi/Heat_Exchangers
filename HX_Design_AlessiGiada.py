import numpy as np
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

# Fins parameters
ef = 0.0005
hf = 0.006
kf = 43
Nf = 12

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

# Outlet Temp. Initial  Guess
T1og = 60 + 273.15
T2og = 40 + 273.15

# Iteration Process
tol = 1e-3
max_iter = 100

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

    # Fin Efficiencies
    m = np.sqrt((2 * a2) / (kf * ef))
    etaf = np.tanh(m * hf) / (m * hf)
    etao = 1 - (Af / Ao) * (1 - etaf)

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


# Print Results
print("=== Results ===")
print(f"Heat Delivered: Q = {Q:.2f} W")
print(f"Water Outlet Temperature: {T1o - 273.15:.2f} °C")
print(f"Air Outlet Temperature: {T2o - 273.15:.2f} °C")
print(f"Water Outlet Velocity: {v1o:.2f} m/s")
print(f"Air Outlet Velocity: {v2o:.2f} m/s")
print(f"Water Pressure Drop: {Dp1:.2f} Pa")
print(f"Air Pressure Drop: {Dp2:.2f} Pa")
print(f"Number of Iterations: {iteration + 1}")
print(f"Water Resistance: {1/a1:.7f}")
print(f"Air Resistance: {1/a2:.7f}")
print(f"Wall Resistance: {ew/kw:.7f}")
print(f"Fin efficiency: {etaf:.4f}")
print(f"Reynolds number: Water = {Re1:.4f}, Air = {Re2:.4f}")
print(f"Prandtl Number: Water = {Pr1:.4f}, Air = {Pr2:.4f}")
print(f"Uo: {Uo:.4f}")
print(f"Mass flow rate of air: m_air = {mf2:.5f}")