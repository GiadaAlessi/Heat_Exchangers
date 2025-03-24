import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from CoolProp.CoolProp import *
import matplotlib.patches as mpatches

#____________________USEFUL FUNCTIONS____________________________________
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from CoolProp.CoolProp import *
import matplotlib.patches as mpatches

#____________________USEFUL FUNCTIONS____________________________________

def toCelsius(T):
    return T-273.15


def toKelvin(T) :
    return T+273.15

def calcul_properties(p,T,Tw,fluid,v,Dh) :
    mu = PropsSI("V","P", p, "T",T,fluid)
    muw = PropsSI("V","P", p, "T",Tw,fluid)
    rho = PropsSI("D", "P", p, "T", T, fluid)
    k = PropsSI("L", "P", p, "T", T, fluid)
    Cp = PropsSI("C", "P", p, "T", T, fluid)

    Re = v * Dh * rho / mu
    Pr = mu * Cp / k
    Gz =  Dh * Re * Pr /  L
    return mu,muw,rho,k,Cp,Re,Pr,Gz

def calcul_alpha_Nu(mu,muw,rho,k,Cp,Re,Pr,Dh,Gz,gas) :
    if Re < 2000:
        if Gz >10 :
            C = 1.86
            m = 1/3
            n = 1/3
            K = ((Dh/L)**(1/3))*(mu/muw)**(0.14)
        else :
            C = 3.66
            m = 0
            n = 0
            K = 1
    elif Re > 2000 and Pr>0.6 and Pr<100 and gas == False:
        C = 0.027
        m = 0.8
        n = 0.33
        K = (mu / muw) ** 0.14

    else :
        C = 0.023
        m = 0.8
        n = 0.4
        K = 1
    Nu = C * (Re ** m) * (Pr ** n) * K
    alpha = Nu * k / Dh
    return alpha,Nu

def calcul_ff(Re,rr) :
    if Re<2000 : f = 16/Re
    elif rr <=0.0001 :
        if Re <3*10e4 : f = 0.079*Re**(-0.25)
        else : f = 0.046*Re**(-0.2)
    elif rr == 0.004 :
        if Re <3*10e4 : f = 0.096*Re**(-0.25)
        else : f = 0.078*Re**(-0.2)
    return f

#_________CREATE A DATAFRAME WITH DATA FROM ATENEA______________________
tabexp = pd.read_csv('datos_CF_titca_formatted.dat', sep='\t+', skipinitialspace=True, comment='#', engine='python')
tabexp.to_csv('output_analysis.csv', index=False)

#____________CODE PROFESSOR________________________________
print ("---------------------------------COUNTERFLOW-----------------------------------------------------------------")
#Plot flows
fig, ax = plt.subplots(1,1, squeeze=False)
fig.suptitle('HX practical activity: Counter flow (our case)')

ax[0,0].plot(tabexp["Time(s)"]/60.0, tabexp["SC-1"],'ro')
ax[0,0].plot(tabexp["Time(s)"]/60.0, tabexp["SC-2"],'b*')
ax[0,0].set_xlabel('Time (min)')
ax[0,0].set_ylabel('Volumetric flow (l/min)')
ax[0,0].grid(True)
ax[0,0].legend(['Hot','Cold'])


fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(12, 6))
fig.suptitle('HX practical activity: T for Counter flow (our case)')
curve1, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-1"], 'ro', label='hot-in')
curve2, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-5"], 'y*', label='hot-out')
curve3, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-10"], 'bo', label='cold-in')
curve4, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-6"], 'g*', label='cold-out')
ax[0,0].set_xlabel('Time (min)')
ax[0,0].set_ylabel('Temperature (ºC)')
# Ajouter des ticks pour l'axe x toutes les 5 minutes
ax[0, 0].set_xticks(range(0, int(max(tabexp["Time(s)"] / 60.0)) + 5, 5))
# Activer la grille (lignes verticales et horizontales)
ax[0, 0].grid(True, which='both', linestyle='--', linewidth=0.5)
#ax[0, 0].legend(['hot-in', 'hot-out', 'cold-in', 'cold-out'], loc='upper left', fontsize=10)
zone1 = ax[0, 0].axvspan(tabexp["Time(s)"][44] / 60.0, tabexp["Time(s)"][81] / 60.0, color='blue', alpha=0.2, label='Mode 1')
# Mode 2 : indices 105 à 146
zone2 = ax[0, 0].axvspan(tabexp["Time(s)"][111] / 60.0, tabexp["Time(s)"][140] / 60.0, color='green', alpha=0.2, label='Mode 2')
# Mode 3 : indices 4 à 198
zone3 = ax[0, 0].axvspan(tabexp["Time(s)"][164] / 60.0, tabexp["Time(s)"][198] / 60.0, color='red', alpha=0.2, label='Mode 3')

# Ajouter la première légende pour les courbes
legend_curves = ax[0, 0].legend(handles=[curve1, curve2, curve3, curve4], loc='upper left', fontsize=10, title='Temperatures')
# Ajouter la deuxième légende pour les zones
legend_modes = ax[0, 0].legend(handles=[mpatches.Patch(color='blue', alpha=0.2, label='Mode 1'),
                                        mpatches.Patch(color='green', alpha=0.2, label='Mode 2'),
                                        mpatches.Patch(color='red', alpha=0.2, label='Mode 3')],
                                loc='upper right', fontsize=10, title='Steady state modes')
# Ajouter les deux légendes à l'axe
ax[0, 0].add_artist(legend_curves)
plt.show()


#__________________INPUT DATA__________________________________
fluid_1 = "Water"
fluid_2 = "Water"
Di = 0.016
Do = 0.018
De1 = 0.026
De2 = 0.028
L = 4
rr = 1e-4
lambda_intermediate = 400 #for conduction between
p1i = 101325
p2i = 101325

#_______________MODES DEFINITION__________________________________________
T1i_1 = np.array(tabexp["ST-1"][44:81])
T1i_2 = np.array(tabexp["ST-1"][111:140])
T1i_3 = np.array(tabexp["ST-1"][164:198])
T1o_1 = np.array(tabexp["ST-5"][44:81])
T1o_2 = np.array(tabexp["ST-5"][111:140])
T1o_3 = np.array(tabexp["ST-5"][164:198])

T2i_1 = np.array(tabexp["ST-10"][44:81])
T2i_2 = np.array(tabexp["ST-10"][111:140])
T2i_3 = np.array(tabexp["ST-10"][164:198])
T2o_1 = np.array(tabexp["ST-6"][44:81])
T2o_2 = np.array(tabexp["ST-6"][111:140])
T2o_3 = np.array(tabexp["ST-6"][164:198])

#flow rate m^3/s
Q1_1 = np.array(tabexp["SC-1"][44:81])/60000
Q1_2 = np.array(tabexp["SC-1"][111:140])/60000
Q1_3 = np.array(tabexp["SC-1"][164:198])/60000

Q2_1 = np.array(tabexp["SC-2"][44:81])/60000
Q2_2 = np.array(tabexp["SC-2"][111:140])/60000
Q2_3 = np.array(tabexp["SC-2"][164:198])/60000

T1i_mean = np.array([np.mean(T1i_1),np.mean(T1i_2),np.mean(T1i_3)])
T2i_mean = np.array([np.mean(T2i_1),np.mean(T2i_2),np.mean(T2i_3)])
T1o_mean = np.array([np.mean(T1o_1),np.mean(T1o_2),np.mean(T1o_3)])
T2o_mean = np.array([np.mean(T2o_1),np.mean(T2o_2),np.mean(T2o_3)])

Q1_mean = np.array([np.mean(Q1_1),np.mean(Q1_2),np.mean(Q1_3)])
Q2_mean = np.array([np.mean(Q2_1),np.mean(Q2_2),np.mean(Q2_3)])

#_______________________TEMPERATURE DIAGRAM STEADY STATE FOR 3 MODES___________________________
label_1i = ['ST-1','ST-2','ST-3','ST-4']
label_1o = ['ST-2','ST-3','ST-4','ST-5']
label_2i = ['ST-7','ST-8','ST-9','ST-10']
label_2o = ['ST-6','ST-7','ST-8','ST-9']
stations_hot = ["ST-1", "ST-2", "ST-3", "ST-4", "ST-5"]
stations_cold = ["ST-6", "ST-7", "ST-8", "ST-9", "ST-10"]

T1_1 = np.zeros(3)
T1_2 = np.zeros(3)
T1_3 = np.zeros(3)
T1_4 = np.zeros(3)
T1_5 = np.zeros(3)

T2_1 = np.zeros(3)
T2_2 = np.zeros(3)
T2_3 = np.zeros(3)
T2_4 = np.zeros(3)
T2_5 = np.zeros(3)

T1_1[0] = np.mean(np.array(tabexp['ST-1'][44:81]))
T1_2[0] = np.mean(np.array(tabexp['ST-2'][44:81]))
T1_3[0] = np.mean(np.array(tabexp['ST-3'][44:81]))
T1_4[0] = np.mean(np.array(tabexp['ST-4'][44:81]))
T1_5[0] = np.mean(np.array(tabexp['ST-5'][44:81]))

T2_1[0] = np.mean(np.array(tabexp['ST-6'][44:81]))
T2_2[0] = np.mean(np.array(tabexp['ST-7'][44:81]))
T2_3[0] = np.mean(np.array(tabexp['ST-8'][44:81]))
T2_4[0] = np.mean(np.array(tabexp['ST-9'][44:81]))
T2_5[0] = np.mean(np.array(tabexp['ST-10'][44:81]))

T1_1[1] = np.mean(np.array(tabexp['ST-1'][111:140]))
T1_2[1] = np.mean(np.array(tabexp['ST-2'][111:140]))
T1_3[1] = np.mean(np.array(tabexp['ST-3'][111:140]))
T1_4[1] = np.mean(np.array(tabexp['ST-4'][111:140]))
T1_5[1] = np.mean(np.array(tabexp['ST-5'][111:140]))

T2_1[1] = np.mean(np.array(tabexp['ST-6'][111:140]))
T2_2[1] = np.mean(np.array(tabexp['ST-7'][111:140]))
T2_3[1] = np.mean(np.array(tabexp['ST-8'][111:140]))
T2_4[1] = np.mean(np.array(tabexp['ST-9'][111:140]))
T2_5[1] = np.mean(np.array(tabexp['ST-10'][111:140]))

T1_1[2] = np.mean(np.array(tabexp['ST-1'][164:198]))
T1_2[2] = np.mean(np.array(tabexp['ST-2'][164:198]))
T1_3[2] = np.mean(np.array(tabexp['ST-3'][164:198]))
T1_4[2] = np.mean(np.array(tabexp['ST-4'][164:198]))
T1_5[2] = np.mean(np.array(tabexp['ST-5'][164:198]))

T2_1[2] = np.mean(np.array(tabexp['ST-6'][164:198]))
T2_2[2] = np.mean(np.array(tabexp['ST-7'][164:198]))
T2_3[2] = np.mean(np.array(tabexp['ST-8'][164:198]))
T2_4[2] = np.mean(np.array(tabexp['ST-9'][164:198]))
T2_5[2] = np.mean(np.array(tabexp['ST-10'][164:198]))

# Création d'une figure avec 3 sous-graphiques côte à côte pour les 3 modes
fig, axs = plt.subplots(1, 3, figsize=(15,5))

# Couleurs de fond pour chaque mode (inspirées du code du professeur)
bg_colors = [(0,0,1,0.2), (0,1,0,0.2), (1,0,0,0.2)]  # bleu, vert, rouge avec alpha

positions = [1,2,3,4,5]
T_amb = [20.1,20.1,20.1,20.1,20.1]

for i in range(3):
    # Définir la couleur de fond du subplot
    axs[i].set_facecolor(bg_colors[i])
    hot_temps = [T1_1[i], T1_2[i], T1_3[i], T1_4[i], T1_5[i]]
    cold_temps = [T2_1[i], T2_2[i], T2_3[i], T2_4[i], T2_5[i]]
   
    # Tracé des températures du fluide chaud pour le mode i
    axs[i].plot(positions, [T1_1[i], T1_2[i], T1_3[i], T1_4[i], T1_5[i]],
                "-", color="red", marker="o", label="hot fluid")
    for idx, (x, y) in enumerate(zip(positions, hot_temps)):
        axs[i].text(x+0.1, y, stations_hot[idx], fontsize=8, color="red")

    # Tracé des températures du fluide froid pour le mode i
    axs[i].plot(positions, [T2_1[i], T2_2[i], T2_3[i], T2_4[i], T2_5[i]],
                "-", color="blue", marker="o", label="cold fluid")
    for idx, (x, y) in enumerate(zip(positions, cold_temps)):
        axs[i].text(x+0.1, y, stations_cold[idx], fontsize=8, color="blue")
   
    # Tracé de la température ambiante
    axs[i].plot(positions, T_amb, "--", color="black", label="Ambient temperature")
   
    axs[i].legend()
    axs[i].set_title("Temperature profile, mode {}".format(i+1))
    axs[i].set_xlabel("Position along the heat exchanger")
    axs[i].set_xticks([])  # On enlève les ticks sur x
    axs[i].set_ylim([16,40])
    if i == 0:
        axs[i].set_ylabel("Temperature [°C]")

plt.tight_layout()
plt.savefig('T_modes_comparison.png')
plt.show()


#____________________________SEMI_ANALYTICAL CALCULATION_________________________________
T1o_analytical = np.zeros(3)
T2o_analytical = np.zeros(3)
Q_analytical   = np.zeros(3)
delta_p1_analytical = np.zeros(3)
delta_p2_analytical = np.zeros(3)
p1o = p1i
p2o = p2i

for i in range(3):  # Pour les 3 modes
    # Conditions d'entrée
    T1i = toKelvin(T1i_mean[i])
    T2i = toKelvin(T2i_mean[i])

    # Hypothèses initiales pour démarrer l'itération
    T1o = T1i - 10
    T2o = T2i + 10

    # Surfaces et diamètres hydrauliques
    A1 = np.pi * Di * L
    S1 = np.pi * Di**2 / 4
    Dh1 = Di  # diam. hydraulique pour un tube circulaire simple = Di
    v1 = Q1_mean[i] / S1

    # Annulaire
    # On considère que le second fluide s'écoule dans l'espace annulaire entre Do et De2
    S2 = np.pi * (De1**2 - Do**2) / 4
    Dh2 = 4 * S2 / (np.pi*(De1+Do))
    v2 = Q2_mean[i] / S2
   
    #Mass flow rate
    m1 = Q1_mean[i]*PropsSI('D', 'T', T1i, 'P', p1i, fluid_1)
    m2 = Q2_mean[i]*PropsSI('D', 'T', T2i, 'P', p2i, fluid_2)

    # Boucle d'itérations pour convergence
    tolerance = 1e-12
    for j in range(5000000):
        T1 = np.mean([T1i, T1o])
        T2 = np.mean([T2i, T2o])
        p1 = np.mean([p1i, p1o])
        p2 = np.mean([p2i, p2o])

        Tw = np.mean([T1,T2])

        # Calcul des propriétés pour chaque fluide
        mu1, muw1, rho1, k1, Cp1, Re1, Pr1, Gz1 = calcul_properties(p1, T1, Tw, fluid_1, v1, Dh1)
        mu2, muw2, rho2, k2, Cp2, Re2, Pr2, Gz2 = calcul_properties(p2, T2, Tw, fluid_2, v2, Dh2)

        # Coefficients de transfert interne et externe
        hi, Nui = calcul_alpha_Nu(mu1, muw1, rho1, k1, Cp1, Re1, Pr1, Dh1, Gz1, False)
        ho, Nuo = calcul_alpha_Nu(mu2, muw2, rho2, k2, Cp2, Re2, Pr2, Dh2, Gz2, False)

        # Pas de rugosité interne ajoutée pour le moment, pas de Rfi et Rfo
        Rconv1 = Do/(hi *Di)
        Rcond  = Do*np.pi*np.log(Do/Di)/(2*np.pi*lambda_intermediate)
        Rconv2 = 1/(ho)
        UA = np.pi*Do*L/(Rconv1 + Rcond + Rconv2)

        # Calcul NTU et epsilon pour contre-courant
        Cmin = min(Cp1*m1, Cp2*m2)
        Cmax = max(Cp1*m1, Cp2*m2)
        Qmax = Cmin * (T1i - T2i)

        NTU = UA / Cmin
        Z   = Cmin / Cmax
        epsilon = (1 - np.exp(-NTU*(1-Z)))/(1-Z*np.exp(-NTU*(1-Z)))

        Q = epsilon * Qmax
        T2o_new = T2i + Q/(m2*Cp2)
        T1o_new = T1i - Q/(m1*Cp1)
       
        if abs(T1o_new - T1o) < tolerance and abs(T2o_new - T2o) < tolerance:
            T1o = T1o_new
            T2o = T2o_new

            # Calcul des pertes de charge
            f1 = calcul_ff(Re1, rr)
            tau1 = f1 * rho1 * v1**2 / 2
            delta_p1 = tau1 * A1 / S1
            p1o = p1i - delta_p1

            f2_i = calcul_ff(Re2, rr)
            tau2_i = f2_i * rho2 * v2**2 / 2
            A2_i = np.pi * Do * L

            f2_e = calcul_ff(Re2, rr)
            tau2_e = f2_e * rho2 * v2**2 / 2
            A2_e = np.pi * De2 * L

            delta_p2 = (tau2_i * A2_i + tau2_e * A2_e)/S2
            p2o = p2i - delta_p2
            break

        # Mise à jour
        T1o = T1o_new
        T2o = T2o_new

        # Calcul des pertes de charge
        f1 = calcul_ff(Re1, rr)
        tau1 = f1 * rho1 * v1**2 / 2
        delta_p1 = tau1 * A1 / S1
        p1o = p1i - delta_p1

        f2_i = calcul_ff(Re2, rr)
        tau2_i = f2_i * rho2 * v2**2 / 2
        A2_i = np.pi * Do * L

        f2_e = calcul_ff(Re2, rr)
        tau2_e = f2_e * rho2 * v2**2 / 2
        A2_e = np.pi * De2 * L

        delta_p2 = (tau2_i * A2_i + tau2_e * A2_e)/S2
        p2o = p2i - delta_p2
       
       
    # Stockage des résultats
    Q_analytical[i] = Q
    T1o_analytical[i] = T1o
    T2o_analytical[i] = T2o
    delta_p1_analytical[i] = delta_p1
    delta_p2_analytical[i] = delta_p2


# Création d'un tableau récapitulatif pour les 3 modes
data = {
    'Mode': [1,2,3],
    'T1i(°C)': T1i_mean,
    'T1o_calc(°C)': T1o_analytical - 273.15,
    'T2i(°C)': T2i_mean,
    'T2o_calc(°C)': T2o_analytical - 273.15,
    'Q(W)': Q_analytical,
    'dp1(kPa)': delta_p1_analytical/1000,
    'dp2(kPa)': delta_p2_analytical/1000
}

df_res = pd.DataFrame(data)
print("_____________Semi-analytical calculation______________")
print(df_res)
#___________________EXPERIMENTAL CALCULATION__________________________
Q_released_exp = np.zeros(3)
Q_absorbed_exp = np.zeros(3)
p1_out = p1i
p2_out = p2i

for j in range(3):
    # Convert average inlet/outlet temperatures to Kelvin
    T_hot_in = toKelvin(T1i_mean[j])
    T_cold_in = toKelvin(T2i_mean[j])
    T_hot_out = toKelvin(T1o_mean[j])
    T_cold_out = toKelvin(T2o_mean[j])
   
    # Compute average properties at mean temperatures and pressures
    Cp_hot = PropsSI("C","T", np.mean([T_hot_in, T_hot_out]), "P", np.mean([p1i,p1_out]), fluid_1)
    Cp_cold = PropsSI("C","T", np.mean([T_cold_in,T_cold_out]), "P", np.mean([p2i, p2_out]), fluid_2)
   
    rho_hot = PropsSI("D","T", T_hot_in, "P", p1i, fluid_1)
    rho_cold = PropsSI("D","T", T_cold_in, "P",p2i, fluid_2)
   
    # Mass flow rates based on density and volumetric flow
    # (Previously named deb1_mean, deb2_mean, now Q1_mean, Q2_mean to ensure consistency)
    m_hot = Q1_mean[j]*rho_hot
    m_cold = Q2_mean[j]*rho_cold
   
    # Calculate experimental heat exchange
    # Q_released: hot fluid releases heat, Q_absorbed: cold fluid absorbs heat
    Q_released = m_hot * Cp_hot * (T_hot_in - T_hot_out)
    Q_absorbed = m_cold * Cp_cold * (T_cold_out - T_cold_in)
   
    Q_loss =   Q_released - Q_absorbed
    efficiency = np.abs(Q_absorbed / Q_released)
   
    Q_released_exp[j] = Q_released
    Q_absorbed_exp[j]=Q_absorbed
   
    # Print results for each mode
    print("_______Experimental Case {}____________".format(j+1))
    print("Q_absorbed_exp = {:.2f} W".format(Q_absorbed))
    print("Q_released_exp = {:.2f} W".format(Q_released))
    print("Q_loss_exp = {:.2f} W".format(Q_loss))
    print("Energy ratio (Q_abs/Q_rel) = {:.2f}".format(efficiency))
    print("_________________________________________")

# Now compare experimental and analytical results in a single table
# We already have Q_analytical from the analytical calculations.
# Let's calculate percentage difference.

difference_percent = np.abs(Q_analytical - Q_released_exp)/Q_analytical * 100
difference_percent_abs = np.abs(Q_analytical - Q_absorbed_exp)/Q_analytical * 100

comparison_data = {
    'Mode': [1, 2, 3],
    'Q_exp(W)': Q_released_exp,
    'Q_analytical(W)': Q_analytical,
    'Difference(%) released': difference_percent,
    'Difference(%) absorbed': difference_percent_abs
    }

df_comparison = pd.DataFrame(comparison_data)
print("______Comparison between Experimental and Analytical Results_________")
print(df_comparison)

#___________________LOSSES DISTRIBTION__________________________
def losses():
    # On se base sur les variables globales déjà définies :
    # T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean, p1i, p2i, fluid_1, fluid_2
    # et le tableau tabexp
   
    label_1i = ['ST-1','ST-2','ST-3','ST-4']
    label_1o = ['ST-2','ST-3','ST-4','ST-5']
    label_2i = ['ST-7','ST-8','ST-9','ST-10']
    label_2o = ['ST-6','ST-7','ST-8','ST-9']

    Qlost_all = np.zeros((4,3))
    Qced_all = np.zeros((4,3))
    Qabs_all = np.zeros((4,3))
    n_BC_all = np.zeros((4,3))

    # On dispose déjà de T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean calculés globalement.
    # Ici, on va supposer que la répartition par portion (j) ne sert qu'à illustrer
    # un calcul local, mais on ne va pas recalculer les moyennes de débit. On utilisera Q1_mean et Q2_mean globaux.
    # Les ranges T1i_range_x, etc. sont lus mais on n'en a pas réellement besoin pour le calcul final,
    # car on a déjà les moyennes globales. Si on veut réellement faire un calcul portion par portion,
    # il faudrait redéfinir la logique. Pour éviter les NaN, on utilisera les moyennes globales déjà existantes.

    for j in range(4):
        # On lit les données portion par portion (on pourrait l'enlever si inutile)
        T1i_range_1 = np.array(tabexp[label_1i[j]][44:81])
        T1i_range_2 = np.array(tabexp[label_1i[j]][111:140])
        T1i_range_3 = np.array(tabexp[label_1i[j]][164:198])
        T1o_range_1 = np.array(tabexp[label_1o[j]][44:81])
        T1o_range_2 = np.array(tabexp[label_1o[j]][111:140])
        T1o_range_3 = np.array(tabexp[label_1o[j]][164:198])

        T2i_range_1 = np.array(tabexp[label_2i[j]][44:81])
        T2i_range_2 = np.array(tabexp[label_2i[j]][111:140])
        T2i_range_3 = np.array(tabexp[label_2i[j]][164:198])
        T2o_range_1 = np.array(tabexp[label_2o[j]][44:81])
        T2o_range_2 = np.array(tabexp[label_2o[j]][111:140])
        T2o_range_3 = np.array(tabexp[label_2o[j]][164:198])

        # On calcule les moyennes par mode pour cette portion
        T1i_portion = np.array([np.mean(T1i_range_1), np.mean(T1i_range_2), np.mean(T1i_range_3)])
        T1o_portion = np.array([np.mean(T1o_range_1), np.mean(T1o_range_2), np.mean(T1o_range_3)])
        T2i_portion = np.array([np.mean(T2i_range_1), np.mean(T2i_range_2), np.mean(T2i_range_3)])
        T2o_portion = np.array([np.mean(T2o_range_1), np.mean(T2o_range_2), np.mean(T2o_range_3)])

        p1o = p1i
        p2o = p2i
        for i in range(3):
            # Convertir en Kelvin
            T1i_K = toKelvin(T1i_portion[i])
            T1o_K = toKelvin(T1o_portion[i])
            T2i_K = toKelvin(T2i_portion[i])
            T2o_K = toKelvin(T2o_portion[i])

            # Utiliser Q1_mean[i] et Q2_mean[i] au lieu de deb1_mean[i], deb2_mean[i]
            m1 = Q1_mean[i] * PropsSI("D","T",T1i_K,"P",p1i,fluid_1)
            m2 = Q2_mean[i] * PropsSI("D","T",T2i_K,"P",p2i,fluid_2)

            Cp1 = PropsSI("C","T",np.mean([T1i_K,T1o_K]),"P",np.mean([p1i,p1o]),fluid_1)
            Cp2 = PropsSI("C","T", np.mean([T2i_K,T2o_K]), "P", np.mean([p2i, p2o]), fluid_2)

            Qced = m1*Cp1*(T1i_K-T1o_K)
            Qabs = m2*Cp2*(T2o_K-T2i_K)

            n_BC_all[j][i] = Qabs/Qced
            Qabs_all[j][i] = Qabs
            Qced_all[j][i] = Qced
            Qlost_all[j][i] = Qced - Qabs

            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("T1i = {:.2f} °C, T2o = {:.2f} °C, T1o = {:.2f} °C, T2i = {:.2f} °C".format(
                T1i_portion[i], T2o_portion[i], T1o_portion[i], T2i_portion[i]))

    for j in range(4):
        for i in range(3):
            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("Qabs = {:.2f} W\nQced = {:.2f} W\nQlost = {:.2f} W\nn_BC = {:.2f}".format(
                Qabs_all[j][i],Qced_all[j][i],Qlost_all[j][i],n_BC_all[j][i]))

    N = 4
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)

    val_1 = Qlost_all.T[0]
    val_2 = Qlost_all.T[1]
    val_3 = Qlost_all.T[2]

    # Pour éviter les erreurs, on peut s'assurer qu'il n'y a pas de NaN
    val_1 = np.nan_to_num(val_1, nan=0.0)
    val_2 = np.nan_to_num(val_2, nan=0.0)
    val_3 = np.nan_to_num(val_3, nan=0.0)

    rects1 = ax.bar(ind, val_1, width, color='blue', alpha = 0.4)
    rects2 = ax.bar(ind + width, val_2, width, color='green', alpha = 0.4)
    rects3 = ax.bar(ind + width * 2, val_3, width, color='red', alpha = 0.4)

    for i in range(3):
        print("-------Mode{}-------".format(i+1))
        print("losses = {:.2f}".format(sum(Qlost_all.T[i])))

    ax.set_ylabel('$\dot{Q}_{losses}$ [W]')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(('Portion 1', 'Portion 2', 'Portion 3', 'Portion 4'))
    ax.legend((rects1[0], rects2[0], rects3[0]), ('Mode 1', 'Mode 2', 'Mode 3'), loc='lower left')

    def autolabel(rects):
        for rect in rects:
            h = rect.get_height()
            # Conversion en int seulement si h n'est pas NaN
            if not np.isnan(h):
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    plt.title('Heat losses along the heat exchanger for different stationary modes')
    plt.savefig("losses")
    plt.show()

# Appel
L = losses()

#________________________________TEST______________________________
def losses():
    # On se base sur les variables globales déjà définies :
    # T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean, p1i, p2i, fluid_1, fluid_2
    # et le tableau tabexp
   
    label_1i = ['ST-1','ST-2','ST-3','ST-4']
    label_1o = ['ST-2','ST-3','ST-4','ST-5']
    label_2i = ['ST-7','ST-8','ST-9','ST-10']
    label_2o = ['ST-6','ST-7','ST-8','ST-9']

    Qlost_all = np.zeros((4,3))
    Qced_all = np.zeros((4,3))
    Qabs_all = np.zeros((4,3))
    n_BC_all = np.zeros((4,3))

    for j in range(4):
        T1i_range_1 = np.array(tabexp[label_1i[j]][44:81])
        T1i_range_2 = np.array(tabexp[label_1i[j]][111:140])
        T1i_range_3 = np.array(tabexp[label_1i[j]][164:198])
        T1o_range_1 = np.array(tabexp[label_1o[j]][44:81])
        T1o_range_2 = np.array(tabexp[label_1o[j]][111:140])
        T1o_range_3 = np.array(tabexp[label_1o[j]][164:198])

        T2i_range_1 = np.array(tabexp[label_2i[j]][44:81])
        T2i_range_2 = np.array(tabexp[label_2i[j]][111:140])
        T2i_range_3 = np.array(tabexp[label_2i[j]][164:198])
        T2o_range_1 = np.array(tabexp[label_2o[j]][44:81])
        T2o_range_2 = np.array(tabexp[label_2o[j]][111:140])
        T2o_range_3 = np.array(tabexp[label_2o[j]][164:198])

        T1i_portion = np.array([np.mean(T1i_range_1), np.mean(T1i_range_2), np.mean(T1i_range_3)])
        T1o_portion = np.array([np.mean(T1o_range_1), np.mean(T1o_range_2), np.mean(T1o_range_3)])
        T2i_portion = np.array([np.mean(T2i_range_1), np.mean(T2i_range_2), np.mean(T2i_range_3)])
        T2o_portion = np.array([np.mean(T2o_range_1), np.mean(T2o_range_2), np.mean(T2o_range_3)])

        p1o = p1i
        p2o = p2i
        for i in range(3):
            T1i_K = toKelvin(T1i_portion[i])
            T1o_K = toKelvin(T1o_portion[i])
            T2i_K = toKelvin(T2i_portion[i])
            T2o_K = toKelvin(T2o_portion[i])

            m1 = Q1_mean[i] * PropsSI("D","T",T1i_K,"P",p1i,fluid_1)
            m2 = Q2_mean[i] * PropsSI("D","T",T2i_K,"P",p2i,fluid_2)

            Cp1 = PropsSI("C","T",np.mean([T1i_K,T1o_K]),"P",np.mean([p1i,p1o]),fluid_1)
            Cp2 = PropsSI("C","T", np.mean([T2i_K,T2o_K]), "P", np.mean([p2i, p2o]), fluid_2)

            Qced = m1*Cp1*(T1i_K-T1o_K)
            Qabs = m2*Cp2*(T2o_K-T2i_K)

            n_BC_all[j][i] = Qabs/Qced
            Qabs_all[j][i] = Qabs
            Qced_all[j][i] = Qced
            Qlost_all[j][i] = Qced - Qabs

            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("T1i = {:.2f} °C, T2o = {:.2f} °C, T1o = {:.2f} °C, T2i = {:.2f} °C".format(
                T1i_portion[i], T2o_portion[i], T1o_portion[i], T2i_portion[i]))

    for j in range(4):
        for i in range(3):
            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("Qabs = {:.2f} W\nQced = {:.2f} W\nQlost = {:.2f} W\nn_BC = {:.2f}".format(
                Qabs_all[j][i],Qced_all[j][i],Qlost_all[j][i],n_BC_all[j][i]))

    return Qabs_all, Qced_all, Qlost_all, n_BC_all

# Appel de la fonction losses()
Qabs_all, Qced_all, Qlost_all, n_BC_all = losses()

# Maintenant que Qced_all, Qabs_all, Qlost_all, n_BC_all sont définis, on peut tracer les graphiques.

# 1) Graphiques Q_loss, Q_absorbed et Q_released pour chaque mode en fonction de la portion
fig, axs = plt.subplots(1, 3, figsize=(15,5))

bg_colors = [(0,0,1,0.2), (0,1,0,0.2), (1,0,0,0.2)]
positions = np.array([1,2,3,4])

for i in range(3):
    axs[i].set_facecolor(bg_colors[i])
   
    axs[i].plot(positions, Qced_all[:,i], 'o-', color='red', label='Q_released (hot)')
    axs[i].plot(positions, Qabs_all[:,i], 's-', color='blue', label='Q_absorbed (cold)')
    axs[i].plot(positions, Qlost_all[:,i], 'd-', color='black', label='Q_loss')
   
    axs[i].set_xlabel("Portion")
    axs[i].set_title("Mode {}".format(i+1))
    axs[i].grid(True, linestyle='--', linewidth=0.5)
    if i == 0:
        axs[i].set_ylabel("Power [W]")
    axs[i].legend()

plt.tight_layout()
plt.savefig('Q_comparison_per_portion.png')
plt.show()

# 2) Graphique Q_loss moyen en fonction du temps
# Les indices des modes sont déjà définis : mode1_indices, mode2_indices, mode3_indices
# On utilise la fonction calc_Q_loss_time déjà définie.

def calc_Q_loss_time(interval_indices, fluid_1, fluid_2):
    Q_loss_time = []
    times = tabexp["Time(s)"][interval_indices]/60.0  # en minutes
    for idx in interval_indices:
        # Températures instantanées
        T1i_inst = toKelvin(tabexp["ST-1"][idx])
        T1o_inst = toKelvin(tabexp["ST-5"][idx])
        T2i_inst = toKelvin(tabexp["ST-10"][idx])
        T2o_inst = toKelvin(tabexp["ST-6"][idx])
       
        # Débits volumiques instantanés
        Q1_inst = tabexp["SC-1"][idx]/60000.0
        Q2_inst = tabexp["SC-2"][idx]/60000.0
       
        # Propriétés
        p1_mean = p1i  # On suppose pression constante
        p2_mean = p2i
       
        T_hot_mean = np.mean([T1i_inst,T1o_inst])
        T_cold_mean = np.mean([T2i_inst,T2o_inst])
       
        Cp_hot = PropsSI("C","T",T_hot_mean,"P",p1_mean,fluid_1)
        Cp_cold = PropsSI("C","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        rho_hot = PropsSI("D","T",T_hot_mean,"P",p1_mean,fluid_1)
        rho_cold = PropsSI("D","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        m_hot = Q1_inst * rho_hot
        m_cold = Q2_inst * rho_cold
       
        Qced_inst = m_hot*Cp_hot*(T1i_inst - T1o_inst)
        Qabs_inst = m_cold*Cp_cold*(T2o_inst - T2i_inst)
       
        Q_loss_inst = Qced_inst - Qabs_inst
        Q_loss_time.append(Q_loss_inst)
   
    return times, np.array(Q_loss_time)


times1, Q_loss_1 = calc_Q_loss_time(range(44,81), fluid_1, fluid_2)
times2, Q_loss_2 = calc_Q_loss_time(range(111,140), fluid_1, fluid_2)
times3, Q_loss_3 = calc_Q_loss_time(range(164,198), fluid_1, fluid_2)

fig, ax = plt.subplots()
ax.plot(times1, Q_loss_1, '-', color='orange', label='Q_loss Mode 1')
ax.plot(times2, Q_loss_2, '-', color='blue', label='Q_loss Mode 2')
ax.plot(times3, Q_loss_3, '-', color='purple', label='Q_loss Mode 3')
ax.set_xlabel('Time [min]')
ax.set_ylabel('Q_loss [W]')
ax.set_title('Q_loss over time for different modes')
ax.grid(True, linestyle='--', linewidth=0.5)
ax.legend()
plt.savefig('Q_loss_over_time.png')
plt.show()


# Recalcul des Q_abs, Q_rel, Q_loss sur toute la durée
def calc_Q_values_time(tabexp, fluid_1, fluid_2, p1i, p2i):
    times = tabexp["Time(s)"].values/60.0  # en minutes
    Q_abs_time = []
    Q_rel_time = []
    Q_loss_time = []
    for idx in range(len(tabexp)):
        T1i_inst = toKelvin(tabexp["ST-1"][idx])
        T1o_inst = toKelvin(tabexp["ST-5"][idx])
        T2i_inst = toKelvin(tabexp["ST-10"][idx])
        T2o_inst = toKelvin(tabexp["ST-6"][idx])
       
        Q1_inst = tabexp["SC-1"][idx]/60000.0
        Q2_inst = tabexp["SC-2"][idx]/60000.0
       
        p1_mean = p1i
        p2_mean = p2i
       
        T_hot_mean = np.mean([T1i_inst,T1o_inst])
        T_cold_mean = np.mean([T2i_inst,T2o_inst])
       
        # Propriétés moyennes
        Cp_hot = PropsSI("C","T",T_hot_mean,"P",p1_mean,fluid_1)
        Cp_cold = PropsSI("C","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        rho_hot = PropsSI("D","T",T_hot_mean,"P",p1_mean,fluid_1)
        rho_cold = PropsSI("D","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        m_hot = Q1_inst * rho_hot
        m_cold = Q2_inst * rho_cold
       
        Q_released_inst = m_hot * Cp_hot * (T1i_inst - T1o_inst)
        Q_absorbed_inst = m_cold * Cp_cold * (T2o_inst - T2i_inst)
        Q_loss_inst = Q_released_inst - Q_absorbed_inst
       
        Q_abs_time.append(Q_absorbed_inst)
        Q_rel_time.append(Q_released_inst)
        Q_loss_time.append(Q_loss_inst)
       
    return times, np.array(Q_abs_time), np.array(Q_rel_time), np.array(Q_loss_time)

times, Q_abs_full, Q_rel_full, Q_loss_full = calc_Q_values_time(tabexp, fluid_1, fluid_2, p1i, p2i)

fig, ax = plt.subplots(figsize=(12,6))
# Tracé des courbes avec nouvelles couleurs
line_abs, = ax.plot(times, Q_abs_full, '-', color='black', label='Q_absorbed')
line_rel, = ax.plot(times, Q_rel_full, '-', color='magenta', label='Q_released')
line_loss,= ax.plot(times, Q_loss_full, '-', color='orange', label='Q_loss')

ax.set_xlabel('Time [min]')
ax.set_ylabel('Power [W]')
ax.set_title('Q_absorbed, Q_released and Q_loss over the entire experiment')
ax.grid(True, linestyle='--', linewidth=0.5)

# Ajout des zones pour les modes, mêmes couleurs (blue, green, red) en alpha=0.2 sans hachures
# Mode 1 : indices 32 à 85
ax.axvspan(times[44], times[81], facecolor='blue', alpha=0.2, label='Mode 1')
# Mode 2 : indices 111 à 143
ax.axvspan(times[111], times[140], facecolor='green', alpha=0.2, label='Mode 2')
# Mode 3 : indices 162 à 198
ax.axvspan(times[164], times[198], facecolor='red', alpha=0.2, label='Mode 3')

# Gestion des légendes
# La première légende pour les courbes Q
legend_curves = ax.legend(handles=[line_abs, line_rel, line_loss], loc='upper left', title='Heat flows')
ax.add_artist(legend_curves)

# Deuxième légende pour les modes
from matplotlib.patches import Patch
legend_modes = [Patch(facecolor='blue', alpha=0.2, label='Mode 1'),
                Patch(facecolor='green', alpha=0.2, label='Mode 2'),
                Patch(facecolor='red', alpha=0.2, label='Mode 3')]
ax.legend(handles=legend_modes, loc='upper right', title='Modes')

# Ajout de petites annotations sur le graphique pour identifier les courbes
ax.text(times[-1]*0.7, np.max(Q_abs_full)*0.9, "Q_absorbed", color='black', fontsize=10)
ax.text(times[-1]*0.7, np.max(Q_rel_full)*0.8, "Q_released", color='magenta', fontsize=10)
ax.text(times[-1]*0.7, np.max(Q_loss_full)*0.7, "Q_loss", color='orange', fontsize=10)

plt.tight_layout()
plt.savefig('Q_all_over_time_with_modes_colored.png')
plt.show()
def toCelsius(T):
    return T-273.15


def toKelvin(T) :
    return T+273.15

def calcul_properties(p,T,Tw,fluid,v,Dh) :
    mu = PropsSI("V","P", p, "T",T,fluid)
    muw = PropsSI("V","P", p, "T",Tw,fluid)
    rho = PropsSI("D", "P", p, "T", T, fluid)
    k = PropsSI("L", "P", p, "T", T, fluid)
    Cp = PropsSI("C", "P", p, "T", T, fluid)

    Re = v * Dh * rho / mu
    Pr = mu * Cp / k
    Gz =  Dh * Re * Pr /  L
    return mu,muw,rho,k,Cp,Re,Pr,Gz

def calcul_alpha_Nu(mu,muw,rho,k,Cp,Re,Pr,Dh,Gz,gas) :
    if Re < 2000:
        if Gz >10 :
            C = 1.86
            m = 1/3
            n = 1/3
            K = ((Dh/L)**(1/3))*(mu/muw)**(0.14)
        else :
            C = 3.66
            m = 0
            n = 0
            K = 1
    elif Re > 2000 and Pr>0.6 and Pr<100 and gas == False:
        C = 0.027
        m = 0.8
        n = 0.33
        K = (mu / muw) ** 0.14

    else :
        C = 0.023
        m = 0.8
        n = 0.4
        K = 1
    Nu = C * (Re ** m) * (Pr ** n) * K
    alpha = Nu * k / Dh
    return alpha,Nu

def calcul_ff(Re,rr) :
    if Re<2000 : f = 16/Re
    elif rr <=0.0001 :
        if Re <3*10e4 : f = 0.079*Re**(-0.25)
        else : f = 0.046*Re**(-0.2)
    elif rr == 0.004 :
        if Re <3*10e4 : f = 0.096*Re**(-0.25)
        else : f = 0.078*Re**(-0.2)
    return f

#_________CREATE A DATAFRAME WITH DATA FROM ATENEA______________________
tabexp = pd.read_csv('datos_CF_titca_formatted.dat', sep='\t+', skipinitialspace=True, comment='#', engine='python')
tabexp.to_csv('output_analysis.csv', index=False)

#____________CODE PROFESSOR________________________________
print ("---------------------------------COUNTERFLOW-----------------------------------------------------------------")
#Plot flows
fig, ax = plt.subplots(1,1, squeeze=False)
fig.suptitle('HX practical activity: Counter flow (our case)')

ax[0,0].plot(tabexp["Time(s)"]/60.0, tabexp["SC-1"],'ro')
ax[0,0].plot(tabexp["Time(s)"]/60.0, tabexp["SC-2"],'b*')
ax[0,0].set_xlabel('Time (min)')
ax[0,0].set_ylabel('Volumetric flow (l/min)')
ax[0,0].grid(True)
ax[0,0].legend(['Hot','Cold'])


fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(12, 6))
fig.suptitle('HX practical activity: T for Counter flow (our case)')
curve1, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-1"], 'ro', label='hot-in')
curve2, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-5"], 'y*', label='hot-out')
curve3, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-10"], 'bo', label='cold-in')
curve4, = ax[0, 0].plot(tabexp["Time(s)"] / 60.0, tabexp["ST-6"], 'g*', label='cold-out')
ax[0,0].set_xlabel('Time (min)')
ax[0,0].set_ylabel('Temperature (ºC)')
# Ajouter des ticks pour l'axe x toutes les 5 minutes
ax[0, 0].set_xticks(range(0, int(max(tabexp["Time(s)"] / 60.0)) + 5, 5))
# Activer la grille (lignes verticales et horizontales)
ax[0, 0].grid(True, which='both', linestyle='--', linewidth=0.5)
#ax[0, 0].legend(['hot-in', 'hot-out', 'cold-in', 'cold-out'], loc='upper left', fontsize=10)
zone1 = ax[0, 0].axvspan(tabexp["Time(s)"][44] / 60.0, tabexp["Time(s)"][81] / 60.0, color='blue', alpha=0.2, label='Mode 1')
# Mode 2 : indices 105 à 146
zone2 = ax[0, 0].axvspan(tabexp["Time(s)"][111] / 60.0, tabexp["Time(s)"][140] / 60.0, color='green', alpha=0.2, label='Mode 2')
# Mode 3 : indices 4 à 198
zone3 = ax[0, 0].axvspan(tabexp["Time(s)"][164] / 60.0, tabexp["Time(s)"][198] / 60.0, color='red', alpha=0.2, label='Mode 3')

# Ajouter la première légende pour les courbes
legend_curves = ax[0, 0].legend(handles=[curve1, curve2, curve3, curve4], loc='upper left', fontsize=10, title='Temperatures')
# Ajouter la deuxième légende pour les zones
legend_modes = ax[0, 0].legend(handles=[mpatches.Patch(color='blue', alpha=0.2, label='Mode 1'),
                                        mpatches.Patch(color='green', alpha=0.2, label='Mode 2'),
                                        mpatches.Patch(color='red', alpha=0.2, label='Mode 3')],
                                loc='upper right', fontsize=10, title='Steady state modes')
# Ajouter les deux légendes à l'axe
ax[0, 0].add_artist(legend_curves)
plt.show()


#__________________INPUT DATA__________________________________
fluid_1 = "Water"
fluid_2 = "Water"
Di = 0.016
Do = 0.018
De1 = 0.026
De2 = 0.028
L = 4
rr = 1e-4
lambda_intermediate = 400 #for conduction between
p1i = 101325
p2i = 101325

#_______________MODES DEFINITION__________________________________________
T1i_1 = np.array(tabexp["ST-1"][44:81])
T1i_2 = np.array(tabexp["ST-1"][111:140])
T1i_3 = np.array(tabexp["ST-1"][164:198])
T1o_1 = np.array(tabexp["ST-5"][44:81])
T1o_2 = np.array(tabexp["ST-5"][111:140])
T1o_3 = np.array(tabexp["ST-5"][164:198])

T2i_1 = np.array(tabexp["ST-10"][44:81])
T2i_2 = np.array(tabexp["ST-10"][111:140])
T2i_3 = np.array(tabexp["ST-10"][164:198])
T2o_1 = np.array(tabexp["ST-6"][44:81])
T2o_2 = np.array(tabexp["ST-6"][111:140])
T2o_3 = np.array(tabexp["ST-6"][164:198])

#flow rate m^3/s
Q1_1 = np.array(tabexp["SC-1"][44:81])/60000
Q1_2 = np.array(tabexp["SC-1"][111:140])/60000
Q1_3 = np.array(tabexp["SC-1"][164:198])/60000

Q2_1 = np.array(tabexp["SC-2"][44:81])/60000
Q2_2 = np.array(tabexp["SC-2"][111:140])/60000
Q2_3 = np.array(tabexp["SC-2"][164:198])/60000

T1i_mean = np.array([np.mean(T1i_1),np.mean(T1i_2),np.mean(T1i_3)])
T2i_mean = np.array([np.mean(T2i_1),np.mean(T2i_2),np.mean(T2i_3)])
T1o_mean = np.array([np.mean(T1o_1),np.mean(T1o_2),np.mean(T1o_3)])
T2o_mean = np.array([np.mean(T2o_1),np.mean(T2o_2),np.mean(T2o_3)])

Q1_mean = np.array([np.mean(Q1_1),np.mean(Q1_2),np.mean(Q1_3)])
Q2_mean = np.array([np.mean(Q2_1),np.mean(Q2_2),np.mean(Q2_3)])

#_______________________TEMPERATURE DIAGRAM STEADY STATE FOR 3 MODES___________________________
label_1i = ['ST-1','ST-2','ST-3','ST-4']
label_1o = ['ST-2','ST-3','ST-4','ST-5']
label_2i = ['ST-7','ST-8','ST-9','ST-10']
label_2o = ['ST-6','ST-7','ST-8','ST-9']
stations_hot = ["ST-1", "ST-2", "ST-3", "ST-4", "ST-5"]
stations_cold = ["ST-6", "ST-7", "ST-8", "ST-9", "ST-10"]

T1_1 = np.zeros(3)
T1_2 = np.zeros(3)
T1_3 = np.zeros(3)
T1_4 = np.zeros(3)
T1_5 = np.zeros(3)

T2_1 = np.zeros(3)
T2_2 = np.zeros(3)
T2_3 = np.zeros(3)
T2_4 = np.zeros(3)
T2_5 = np.zeros(3)

T1_1[0] = np.mean(np.array(tabexp['ST-1'][44:81]))
T1_2[0] = np.mean(np.array(tabexp['ST-2'][44:81]))
T1_3[0] = np.mean(np.array(tabexp['ST-3'][44:81]))
T1_4[0] = np.mean(np.array(tabexp['ST-4'][44:81]))
T1_5[0] = np.mean(np.array(tabexp['ST-5'][44:81]))

T2_1[0] = np.mean(np.array(tabexp['ST-6'][44:81]))
T2_2[0] = np.mean(np.array(tabexp['ST-7'][44:81]))
T2_3[0] = np.mean(np.array(tabexp['ST-8'][44:81]))
T2_4[0] = np.mean(np.array(tabexp['ST-9'][44:81]))
T2_5[0] = np.mean(np.array(tabexp['ST-10'][44:81]))

T1_1[1] = np.mean(np.array(tabexp['ST-1'][111:140]))
T1_2[1] = np.mean(np.array(tabexp['ST-2'][111:140]))
T1_3[1] = np.mean(np.array(tabexp['ST-3'][111:140]))
T1_4[1] = np.mean(np.array(tabexp['ST-4'][111:140]))
T1_5[1] = np.mean(np.array(tabexp['ST-5'][111:140]))

T2_1[1] = np.mean(np.array(tabexp['ST-6'][111:140]))
T2_2[1] = np.mean(np.array(tabexp['ST-7'][111:140]))
T2_3[1] = np.mean(np.array(tabexp['ST-8'][111:140]))
T2_4[1] = np.mean(np.array(tabexp['ST-9'][111:140]))
T2_5[1] = np.mean(np.array(tabexp['ST-10'][111:140]))

T1_1[2] = np.mean(np.array(tabexp['ST-1'][164:198]))
T1_2[2] = np.mean(np.array(tabexp['ST-2'][164:198]))
T1_3[2] = np.mean(np.array(tabexp['ST-3'][164:198]))
T1_4[2] = np.mean(np.array(tabexp['ST-4'][164:198]))
T1_5[2] = np.mean(np.array(tabexp['ST-5'][164:198]))

T2_1[2] = np.mean(np.array(tabexp['ST-6'][164:198]))
T2_2[2] = np.mean(np.array(tabexp['ST-7'][164:198]))
T2_3[2] = np.mean(np.array(tabexp['ST-8'][164:198]))
T2_4[2] = np.mean(np.array(tabexp['ST-9'][164:198]))
T2_5[2] = np.mean(np.array(tabexp['ST-10'][164:198]))

# Création d'une figure avec 3 sous-graphiques côte à côte pour les 3 modes
fig, axs = plt.subplots(1, 3, figsize=(15,5))

# Couleurs de fond pour chaque mode (inspirées du code du professeur)
bg_colors = [(0,0,1,0.2), (0,1,0,0.2), (1,0,0,0.2)]  # bleu, vert, rouge avec alpha

positions = [1,2,3,4,5]
T_amb = [20.1,20.1,20.1,20.1,20.1]

for i in range(3):
    # Définir la couleur de fond du subplot
    axs[i].set_facecolor(bg_colors[i])
    hot_temps = [T1_1[i], T1_2[i], T1_3[i], T1_4[i], T1_5[i]]
    cold_temps = [T2_1[i], T2_2[i], T2_3[i], T2_4[i], T2_5[i]]
   
    # Tracé des températures du fluide chaud pour le mode i
    axs[i].plot(positions, [T1_1[i], T1_2[i], T1_3[i], T1_4[i], T1_5[i]],
                "-", color="red", marker="o", label="hot fluid")
    for idx, (x, y) in enumerate(zip(positions, hot_temps)):
        axs[i].text(x+0.1, y, stations_hot[idx], fontsize=8, color="red")

    # Tracé des températures du fluide froid pour le mode i
    axs[i].plot(positions, [T2_1[i], T2_2[i], T2_3[i], T2_4[i], T2_5[i]],
                "-", color="blue", marker="o", label="cold fluid")
    for idx, (x, y) in enumerate(zip(positions, cold_temps)):
        axs[i].text(x+0.1, y, stations_cold[idx], fontsize=8, color="blue")
   
    # Tracé de la température ambiante
    axs[i].plot(positions, T_amb, "--", color="black", label="Ambient temperature")
   
    axs[i].legend()
    axs[i].set_title("Temperature profile, mode {}".format(i+1))
    axs[i].set_xlabel("Position along the heat exchanger")
    axs[i].set_xticks([])  # On enlève les ticks sur x
    axs[i].set_ylim([16,40])
    if i == 0:
        axs[i].set_ylabel("Temperature [°C]")

plt.tight_layout()
plt.savefig('T_modes_comparison.png')
plt.show()


#____________________________SEMI_ANALYTICAL CALCULATION_________________________________
T1o_analytical = np.zeros(3)
T2o_analytical = np.zeros(3)
Q_analytical   = np.zeros(3)
delta_p1_analytical = np.zeros(3)
delta_p2_analytical = np.zeros(3)
p1o = p1i
p2o = p2i

for i in range(3):  # Pour les 3 modes
    # Conditions d'entrée
    T1i = toKelvin(T1i_mean[i])
    T2i = toKelvin(T2i_mean[i])

    # Hypothèses initiales pour démarrer l'itération
    T1o = T1i - 10
    T2o = T2i + 10

    # Surfaces et diamètres hydrauliques
    A1 = np.pi * Di * L
    S1 = np.pi * Di**2 / 4
    Dh1 = Di  # diam. hydraulique pour un tube circulaire simple = Di
    v1 = Q1_mean[i] / S1

    # Annulaire
    # On considère que le second fluide s'écoule dans l'espace annulaire entre Do et De2
    S2 = np.pi * (De1**2 - Do**2) / 4
    Dh2 = 4 * S2 / (np.pi*(De1+Do))
    v2 = Q2_mean[i] / S2
   
    #Mass flow rate
    m1 = Q1_mean[i]*PropsSI('D', 'T', T1i, 'P', p1i, fluid_1)
    m2 = Q2_mean[i]*PropsSI('D', 'T', T2i, 'P', p2i, fluid_2)

    # Boucle d'itérations pour convergence
    tolerance = 1e-12
    for j in range(5000000):
        T1 = np.mean([T1i, T1o])
        T2 = np.mean([T2i, T2o])
        p1 = np.mean([p1i, p1o])
        p2 = np.mean([p2i, p2o])

        Tw = np.mean([T1,T2])

        # Calcul des propriétés pour chaque fluide
        mu1, muw1, rho1, k1, Cp1, Re1, Pr1, Gz1 = calcul_properties(p1, T1, Tw, fluid_1, v1, Dh1)
        mu2, muw2, rho2, k2, Cp2, Re2, Pr2, Gz2 = calcul_properties(p2, T2, Tw, fluid_2, v2, Dh2)

        # Coefficients de transfert interne et externe
        hi, Nui = calcul_alpha_Nu(mu1, muw1, rho1, k1, Cp1, Re1, Pr1, Dh1, Gz1, False)
        ho, Nuo = calcul_alpha_Nu(mu2, muw2, rho2, k2, Cp2, Re2, Pr2, Dh2, Gz2, False)

        # Pas de rugosité interne ajoutée pour le moment, pas de Rfi et Rfo
        Rconv1 = Do/(hi *Di)
        Rcond  = Do*np.pi*np.log(Do/Di)/(2*np.pi*lambda_intermediate)
        Rconv2 = 1/(ho)
        UA = np.pi*Do*L/(Rconv1 + Rcond + Rconv2)

        # Calcul NTU et epsilon pour contre-courant
        Cmin = min(Cp1*m1, Cp2*m2)
        Cmax = max(Cp1*m1, Cp2*m2)
        Qmax = Cmin * (T1i - T2i)

        NTU = UA / Cmin
        Z   = Cmin / Cmax
        epsilon = (1 - np.exp(-NTU*(1-Z)))/(1-Z*np.exp(-NTU*(1-Z)))

        Q = epsilon * Qmax
        T2o_new = T2i + Q/(m2*Cp2)
        T1o_new = T1i - Q/(m1*Cp1)
       
        if abs(T1o_new - T1o) < tolerance and abs(T2o_new - T2o) < tolerance:
            T1o = T1o_new
            T2o = T2o_new

            # Calcul des pertes de charge
            f1 = calcul_ff(Re1, rr)
            tau1 = f1 * rho1 * v1**2 / 2
            delta_p1 = tau1 * A1 / S1
            p1o = p1i - delta_p1

            f2_i = calcul_ff(Re2, rr)
            tau2_i = f2_i * rho2 * v2**2 / 2
            A2_i = np.pi * Do * L

            f2_e = calcul_ff(Re2, rr)
            tau2_e = f2_e * rho2 * v2**2 / 2
            A2_e = np.pi * De2 * L

            delta_p2 = (tau2_i * A2_i + tau2_e * A2_e)/S2
            p2o = p2i - delta_p2
            break

        # Mise à jour
        T1o = T1o_new
        T2o = T2o_new

        # Calcul des pertes de charge
        f1 = calcul_ff(Re1, rr)
        tau1 = f1 * rho1 * v1**2 / 2
        delta_p1 = tau1 * A1 / S1
        p1o = p1i - delta_p1

        f2_i = calcul_ff(Re2, rr)
        tau2_i = f2_i * rho2 * v2**2 / 2
        A2_i = np.pi * Do * L

        f2_e = calcul_ff(Re2, rr)
        tau2_e = f2_e * rho2 * v2**2 / 2
        A2_e = np.pi * De2 * L

        delta_p2 = (tau2_i * A2_i + tau2_e * A2_e)/S2
        p2o = p2i - delta_p2
       
       
    # Stockage des résultats
    Q_analytical[i] = Q
    T1o_analytical[i] = T1o
    T2o_analytical[i] = T2o
    delta_p1_analytical[i] = delta_p1
    delta_p2_analytical[i] = delta_p2


# Création d'un tableau récapitulatif pour les 3 modes
data = {
    'Mode': [1,2,3],
    'T1i(°C)': T1i_mean,
    'T1o_calc(°C)': T1o_analytical - 273.15,
    'T2i(°C)': T2i_mean,
    'T2o_calc(°C)': T2o_analytical - 273.15,
    'Q(W)': Q_analytical,
    'dp1(kPa)': delta_p1_analytical/1000,
    'dp2(kPa)': delta_p2_analytical/1000
}

df_res = pd.DataFrame(data)
print("_____________Semi-analytical calculation______________")
print(df_res)
#___________________EXPERIMENTAL CALCULATION__________________________
Q_released_exp = np.zeros(3)
Q_absorbed_exp = np.zeros(3)
p1_out = p1i
p2_out = p2i

for j in range(3):
    # Convert average inlet/outlet temperatures to Kelvin
    T_hot_in = toKelvin(T1i_mean[j])
    T_cold_in = toKelvin(T2i_mean[j])
    T_hot_out = toKelvin(T1o_mean[j])
    T_cold_out = toKelvin(T2o_mean[j])
   
    # Compute average properties at mean temperatures and pressures
    Cp_hot = PropsSI("C","T", np.mean([T_hot_in, T_hot_out]), "P", np.mean([p1i,p1_out]), fluid_1)
    Cp_cold = PropsSI("C","T", np.mean([T_cold_in,T_cold_out]), "P", np.mean([p2i, p2_out]), fluid_2)
   
    rho_hot = PropsSI("D","T", T_hot_in, "P", p1i, fluid_1)
    rho_cold = PropsSI("D","T", T_cold_in, "P",p2i, fluid_2)
   
    # Mass flow rates based on density and volumetric flow
    # (Previously named deb1_mean, deb2_mean, now Q1_mean, Q2_mean to ensure consistency)
    m_hot = Q1_mean[j]*rho_hot
    m_cold = Q2_mean[j]*rho_cold
   
    # Calculate experimental heat exchange
    # Q_released: hot fluid releases heat, Q_absorbed: cold fluid absorbs heat
    Q_released = m_hot * Cp_hot * (T_hot_in - T_hot_out)
    Q_absorbed = m_cold * Cp_cold * (T_cold_out - T_cold_in)
   
    Q_loss =   Q_released - Q_absorbed
    efficiency = np.abs(Q_absorbed / Q_released)
   
    Q_released_exp[j] = Q_released
    Q_absorbed_exp[j]=Q_absorbed
   
    # Print results for each mode
    print("_______Experimental Case {}____________".format(j+1))
    print("Q_absorbed_exp = {:.2f} W".format(Q_absorbed))
    print("Q_released_exp = {:.2f} W".format(Q_released))
    print("Q_loss_exp = {:.2f} W".format(Q_loss))
    print("Energy ratio (Q_abs/Q_rel) = {:.2f}".format(efficiency))
    print("_________________________________________")

# Now compare experimental and analytical results in a single table
# We already have Q_analytical from the analytical calculations.
# Let's calculate percentage difference.

difference_percent = np.abs(Q_analytical - Q_released_exp)/Q_analytical * 100
difference_percent_abs = np.abs(Q_analytical - Q_absorbed_exp)/Q_analytical * 100

comparison_data = {
    'Mode': [1, 2, 3],
    'Q_exp(W)': Q_released_exp,
    'Q_analytical(W)': Q_analytical,
    'Difference(%) released': difference_percent,
    'Difference(%) absorbed': difference_percent_abs
    }

df_comparison = pd.DataFrame(comparison_data)
print("______Comparison between Experimental and Analytical Results_________")
print(df_comparison)

#___________________LOSSES DISTRIBTION__________________________
def losses():
    # On se base sur les variables globales déjà définies :
    # T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean, p1i, p2i, fluid_1, fluid_2
    # et le tableau tabexp
   
    label_1i = ['ST-1','ST-2','ST-3','ST-4']
    label_1o = ['ST-2','ST-3','ST-4','ST-5']
    label_2i = ['ST-7','ST-8','ST-9','ST-10']
    label_2o = ['ST-6','ST-7','ST-8','ST-9']

    Qlost_all = np.zeros((4,3))
    Qced_all = np.zeros((4,3))
    Qabs_all = np.zeros((4,3))
    n_BC_all = np.zeros((4,3))

    # On dispose déjà de T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean calculés globalement.
    # Ici, on va supposer que la répartition par portion (j) ne sert qu'à illustrer
    # un calcul local, mais on ne va pas recalculer les moyennes de débit. On utilisera Q1_mean et Q2_mean globaux.
    # Les ranges T1i_range_x, etc. sont lus mais on n'en a pas réellement besoin pour le calcul final,
    # car on a déjà les moyennes globales. Si on veut réellement faire un calcul portion par portion,
    # il faudrait redéfinir la logique. Pour éviter les NaN, on utilisera les moyennes globales déjà existantes.

    for j in range(4):
        # On lit les données portion par portion (on pourrait l'enlever si inutile)
        T1i_range_1 = np.array(tabexp[label_1i[j]][44:81])
        T1i_range_2 = np.array(tabexp[label_1i[j]][111:140])
        T1i_range_3 = np.array(tabexp[label_1i[j]][164:198])
        T1o_range_1 = np.array(tabexp[label_1o[j]][44:81])
        T1o_range_2 = np.array(tabexp[label_1o[j]][111:140])
        T1o_range_3 = np.array(tabexp[label_1o[j]][164:198])

        T2i_range_1 = np.array(tabexp[label_2i[j]][44:81])
        T2i_range_2 = np.array(tabexp[label_2i[j]][111:140])
        T2i_range_3 = np.array(tabexp[label_2i[j]][164:198])
        T2o_range_1 = np.array(tabexp[label_2o[j]][44:81])
        T2o_range_2 = np.array(tabexp[label_2o[j]][111:140])
        T2o_range_3 = np.array(tabexp[label_2o[j]][164:198])

        # On calcule les moyennes par mode pour cette portion
        T1i_portion = np.array([np.mean(T1i_range_1), np.mean(T1i_range_2), np.mean(T1i_range_3)])
        T1o_portion = np.array([np.mean(T1o_range_1), np.mean(T1o_range_2), np.mean(T1o_range_3)])
        T2i_portion = np.array([np.mean(T2i_range_1), np.mean(T2i_range_2), np.mean(T2i_range_3)])
        T2o_portion = np.array([np.mean(T2o_range_1), np.mean(T2o_range_2), np.mean(T2o_range_3)])

        p1o = p1i
        p2o = p2i
        for i in range(3):
            # Convertir en Kelvin
            T1i_K = toKelvin(T1i_portion[i])
            T1o_K = toKelvin(T1o_portion[i])
            T2i_K = toKelvin(T2i_portion[i])
            T2o_K = toKelvin(T2o_portion[i])

            # Utiliser Q1_mean[i] et Q2_mean[i] au lieu de deb1_mean[i], deb2_mean[i]
            m1 = Q1_mean[i] * PropsSI("D","T",T1i_K,"P",p1i,fluid_1)
            m2 = Q2_mean[i] * PropsSI("D","T",T2i_K,"P",p2i,fluid_2)

            Cp1 = PropsSI("C","T",np.mean([T1i_K,T1o_K]),"P",np.mean([p1i,p1o]),fluid_1)
            Cp2 = PropsSI("C","T", np.mean([T2i_K,T2o_K]), "P", np.mean([p2i, p2o]), fluid_2)

            Qced = m1*Cp1*(T1i_K-T1o_K)
            Qabs = m2*Cp2*(T2o_K-T2i_K)

            n_BC_all[j][i] = Qabs/Qced
            Qabs_all[j][i] = Qabs
            Qced_all[j][i] = Qced
            Qlost_all[j][i] = Qced - Qabs

            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("T1i = {:.2f} °C, T2o = {:.2f} °C, T1o = {:.2f} °C, T2i = {:.2f} °C".format(
                T1i_portion[i], T2o_portion[i], T1o_portion[i], T2i_portion[i]))

    for j in range(4):
        for i in range(3):
            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("Qabs = {:.2f} W\nQced = {:.2f} W\nQlost = {:.2f} W\nn_BC = {:.2f}".format(
                Qabs_all[j][i],Qced_all[j][i],Qlost_all[j][i],n_BC_all[j][i]))

    N = 4
    ind = np.arange(N)  # the x locations for the groups
    width = 0.27       # the width of the bars

    fig = plt.figure()
    ax = fig.add_subplot(111)

    val_1 = Qlost_all.T[0]
    val_2 = Qlost_all.T[1]
    val_3 = Qlost_all.T[2]

    # Pour éviter les erreurs, on peut s'assurer qu'il n'y a pas de NaN
    val_1 = np.nan_to_num(val_1, nan=0.0)
    val_2 = np.nan_to_num(val_2, nan=0.0)
    val_3 = np.nan_to_num(val_3, nan=0.0)

    rects1 = ax.bar(ind, val_1, width, color='blue', alpha = 0.4)
    rects2 = ax.bar(ind + width, val_2, width, color='green', alpha = 0.4)
    rects3 = ax.bar(ind + width * 2, val_3, width, color='red', alpha = 0.4)

    for i in range(3):
        print("-------Mode{}-------".format(i+1))
        print("losses = {:.2f}".format(sum(Qlost_all.T[i])))

    ax.set_ylabel('$\dot{Q}_{losses}$ [W]')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(('Portion 1', 'Portion 2', 'Portion 3', 'Portion 4'))
    ax.legend((rects1[0], rects2[0], rects3[0]), ('Mode 1', 'Mode 2', 'Mode 3'), loc='lower left')

    def autolabel(rects):
        for rect in rects:
            h = rect.get_height()
            # Conversion en int seulement si h n'est pas NaN
            if not np.isnan(h):
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    plt.title('Heat losses along the heat exchanger for different stationary modes')
    plt.savefig("losses")
    plt.show()

# Appel
L = losses()

#________________________________TEST______________________________
def losses():
    # On se base sur les variables globales déjà définies :
    # T1i_mean, T1o_mean, T2i_mean, T2o_mean, Q1_mean, Q2_mean, p1i, p2i, fluid_1, fluid_2
    # et le tableau tabexp
   
    label_1i = ['ST-1','ST-2','ST-3','ST-4']
    label_1o = ['ST-2','ST-3','ST-4','ST-5']
    label_2i = ['ST-7','ST-8','ST-9','ST-10']
    label_2o = ['ST-6','ST-7','ST-8','ST-9']

    Qlost_all = np.zeros((4,3))
    Qced_all = np.zeros((4,3))
    Qabs_all = np.zeros((4,3))
    n_BC_all = np.zeros((4,3))

    for j in range(4):
        T1i_range_1 = np.array(tabexp[label_1i[j]][44:81])
        T1i_range_2 = np.array(tabexp[label_1i[j]][111:140])
        T1i_range_3 = np.array(tabexp[label_1i[j]][164:198])
        T1o_range_1 = np.array(tabexp[label_1o[j]][44:81])
        T1o_range_2 = np.array(tabexp[label_1o[j]][111:140])
        T1o_range_3 = np.array(tabexp[label_1o[j]][164:198])

        T2i_range_1 = np.array(tabexp[label_2i[j]][44:81])
        T2i_range_2 = np.array(tabexp[label_2i[j]][111:140])
        T2i_range_3 = np.array(tabexp[label_2i[j]][164:198])
        T2o_range_1 = np.array(tabexp[label_2o[j]][44:81])
        T2o_range_2 = np.array(tabexp[label_2o[j]][111:140])
        T2o_range_3 = np.array(tabexp[label_2o[j]][164:198])

        T1i_portion = np.array([np.mean(T1i_range_1), np.mean(T1i_range_2), np.mean(T1i_range_3)])
        T1o_portion = np.array([np.mean(T1o_range_1), np.mean(T1o_range_2), np.mean(T1o_range_3)])
        T2i_portion = np.array([np.mean(T2i_range_1), np.mean(T2i_range_2), np.mean(T2i_range_3)])
        T2o_portion = np.array([np.mean(T2o_range_1), np.mean(T2o_range_2), np.mean(T2o_range_3)])

        p1o = p1i
        p2o = p2i
        for i in range(3):
            T1i_K = toKelvin(T1i_portion[i])
            T1o_K = toKelvin(T1o_portion[i])
            T2i_K = toKelvin(T2i_portion[i])
            T2o_K = toKelvin(T2o_portion[i])

            m1 = Q1_mean[i] * PropsSI("D","T",T1i_K,"P",p1i,fluid_1)
            m2 = Q2_mean[i] * PropsSI("D","T",T2i_K,"P",p2i,fluid_2)

            Cp1 = PropsSI("C","T",np.mean([T1i_K,T1o_K]),"P",np.mean([p1i,p1o]),fluid_1)
            Cp2 = PropsSI("C","T", np.mean([T2i_K,T2o_K]), "P", np.mean([p2i, p2o]), fluid_2)

            Qced = m1*Cp1*(T1i_K-T1o_K)
            Qabs = m2*Cp2*(T2o_K-T2i_K)

            n_BC_all[j][i] = Qabs/Qced
            Qabs_all[j][i] = Qabs
            Qced_all[j][i] = Qced
            Qlost_all[j][i] = Qced - Qabs

            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("T1i = {:.2f} °C, T2o = {:.2f} °C, T1o = {:.2f} °C, T2i = {:.2f} °C".format(
                T1i_portion[i], T2o_portion[i], T1o_portion[i], T2i_portion[i]))

    for j in range(4):
        for i in range(3):
            print("=============Portion {}, Mode {}=============".format(j+1,i+1))
            print("Qabs = {:.2f} W\nQced = {:.2f} W\nQlost = {:.2f} W\nn_BC = {:.2f}".format(
                Qabs_all[j][i],Qced_all[j][i],Qlost_all[j][i],n_BC_all[j][i]))

    return Qabs_all, Qced_all, Qlost_all, n_BC_all

# Appel de la fonction losses()
Qabs_all, Qced_all, Qlost_all, n_BC_all = losses()

# Maintenant que Qced_all, Qabs_all, Qlost_all, n_BC_all sont définis, on peut tracer les graphiques.

# 1) Graphiques Q_loss, Q_absorbed et Q_released pour chaque mode en fonction de la portion
fig, axs = plt.subplots(1, 3, figsize=(15,5))

bg_colors = [(0,0,1,0.2), (0,1,0,0.2), (1,0,0,0.2)]
positions = np.array([1,2,3,4])

for i in range(3):
    axs[i].set_facecolor(bg_colors[i])
   
    axs[i].plot(positions, Qced_all[:,i], 'o-', color='red', label='Q_released (hot)')
    axs[i].plot(positions, Qabs_all[:,i], 's-', color='blue', label='Q_absorbed (cold)')
    axs[i].plot(positions, Qlost_all[:,i], 'd-', color='black', label='Q_loss')
   
    axs[i].set_xlabel("Portion")
    axs[i].set_title("Mode {}".format(i+1))
    axs[i].grid(True, linestyle='--', linewidth=0.5)
    if i == 0:
        axs[i].set_ylabel("Power [W]")
    axs[i].legend()

plt.tight_layout()
plt.savefig('Q_comparison_per_portion.png')
plt.show()

# 2) Graphique Q_loss moyen en fonction du temps
# Les indices des modes sont déjà définis : mode1_indices, mode2_indices, mode3_indices
# On utilise la fonction calc_Q_loss_time déjà définie.

def calc_Q_loss_time(interval_indices, fluid_1, fluid_2):
    Q_loss_time = []
    times = tabexp["Time(s)"][interval_indices]/60.0  # en minutes
    for idx in interval_indices:
        # Températures instantanées
        T1i_inst = toKelvin(tabexp["ST-1"][idx])
        T1o_inst = toKelvin(tabexp["ST-5"][idx])
        T2i_inst = toKelvin(tabexp["ST-10"][idx])
        T2o_inst = toKelvin(tabexp["ST-6"][idx])
       
        # Débits volumiques instantanés
        Q1_inst = tabexp["SC-1"][idx]/60000.0
        Q2_inst = tabexp["SC-2"][idx]/60000.0
       
        # Propriétés
        p1_mean = p1i  # On suppose pression constante
        p2_mean = p2i
       
        T_hot_mean = np.mean([T1i_inst,T1o_inst])
        T_cold_mean = np.mean([T2i_inst,T2o_inst])
       
        Cp_hot = PropsSI("C","T",T_hot_mean,"P",p1_mean,fluid_1)
        Cp_cold = PropsSI("C","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        rho_hot = PropsSI("D","T",T_hot_mean,"P",p1_mean,fluid_1)
        rho_cold = PropsSI("D","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        m_hot = Q1_inst * rho_hot
        m_cold = Q2_inst * rho_cold
       
        Qced_inst = m_hot*Cp_hot*(T1i_inst - T1o_inst)
        Qabs_inst = m_cold*Cp_cold*(T2o_inst - T2i_inst)
       
        Q_loss_inst = Qced_inst - Qabs_inst
        Q_loss_time.append(Q_loss_inst)
   
    return times, np.array(Q_loss_time)


times1, Q_loss_1 = calc_Q_loss_time(range(44,81), fluid_1, fluid_2)
times2, Q_loss_2 = calc_Q_loss_time(range(111,140), fluid_1, fluid_2)
times3, Q_loss_3 = calc_Q_loss_time(range(164,198), fluid_1, fluid_2)

fig, ax = plt.subplots()
ax.plot(times1, Q_loss_1, '-', color='orange', label='Q_loss Mode 1')
ax.plot(times2, Q_loss_2, '-', color='blue', label='Q_loss Mode 2')
ax.plot(times3, Q_loss_3, '-', color='purple', label='Q_loss Mode 3')
ax.set_xlabel('Time [min]')
ax.set_ylabel('Q_loss [W]')
ax.set_title('Q_loss over time for different modes')
ax.grid(True, linestyle='--', linewidth=0.5)
ax.legend()
plt.savefig('Q_loss_over_time.png')
plt.show()


# Recalcul des Q_abs, Q_rel, Q_loss sur toute la durée
def calc_Q_values_time(tabexp, fluid_1, fluid_2, p1i, p2i):
    times = tabexp["Time(s)"].values/60.0  # en minutes
    Q_abs_time = []
    Q_rel_time = []
    Q_loss_time = []
    for idx in range(len(tabexp)):
        T1i_inst = toKelvin(tabexp["ST-1"][idx])
        T1o_inst = toKelvin(tabexp["ST-5"][idx])
        T2i_inst = toKelvin(tabexp["ST-10"][idx])
        T2o_inst = toKelvin(tabexp["ST-6"][idx])
       
        Q1_inst = tabexp["SC-1"][idx]/60000.0
        Q2_inst = tabexp["SC-2"][idx]/60000.0
       
        p1_mean = p1i
        p2_mean = p2i
       
        T_hot_mean = np.mean([T1i_inst,T1o_inst])
        T_cold_mean = np.mean([T2i_inst,T2o_inst])
       
        # Propriétés moyennes
        Cp_hot = PropsSI("C","T",T_hot_mean,"P",p1_mean,fluid_1)
        Cp_cold = PropsSI("C","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        rho_hot = PropsSI("D","T",T_hot_mean,"P",p1_mean,fluid_1)
        rho_cold = PropsSI("D","T",T_cold_mean,"P",p2_mean,fluid_2)
       
        m_hot = Q1_inst * rho_hot
        m_cold = Q2_inst * rho_cold
       
        Q_released_inst = m_hot * Cp_hot * (T1i_inst - T1o_inst)
        Q_absorbed_inst = m_cold * Cp_cold * (T2o_inst - T2i_inst)
        Q_loss_inst = Q_released_inst - Q_absorbed_inst
       
        Q_abs_time.append(Q_absorbed_inst)
        Q_rel_time.append(Q_released_inst)
        Q_loss_time.append(Q_loss_inst)
       
    return times, np.array(Q_abs_time), np.array(Q_rel_time), np.array(Q_loss_time)

times, Q_abs_full, Q_rel_full, Q_loss_full = calc_Q_values_time(tabexp, fluid_1, fluid_2, p1i, p2i)

fig, ax = plt.subplots(figsize=(12,6))
# Tracé des courbes avec nouvelles couleurs
line_abs, = ax.plot(times, Q_abs_full, '-', color='black', label='Q_absorbed')
line_rel, = ax.plot(times, Q_rel_full, '-', color='magenta', label='Q_released')
line_loss,= ax.plot(times, Q_loss_full, '-', color='orange', label='Q_loss')

ax.set_xlabel('Time [min]')
ax.set_ylabel('Power [W]')
ax.set_title('Q_absorbed, Q_released and Q_loss over the entire experiment')
ax.grid(True, linestyle='--', linewidth=0.5)

# Ajout des zones pour les modes, mêmes couleurs (blue, green, red) en alpha=0.2 sans hachures
# Mode 1 : indices 32 à 85
ax.axvspan(times[44], times[81], facecolor='blue', alpha=0.2, label='Mode 1')
# Mode 2 : indices 111 à 143
ax.axvspan(times[111], times[140], facecolor='green', alpha=0.2, label='Mode 2')
# Mode 3 : indices 162 à 198
ax.axvspan(times[164], times[198], facecolor='red', alpha=0.2, label='Mode 3')

# Gestion des légendes
# La première légende pour les courbes Q
legend_curves = ax.legend(handles=[line_abs, line_rel, line_loss], loc='upper left', title='Heat flows')
ax.add_artist(legend_curves)

# Deuxième légende pour les modes
from matplotlib.patches import Patch
legend_modes = [Patch(facecolor='blue', alpha=0.2, label='Mode 1'),
                Patch(facecolor='green', alpha=0.2, label='Mode 2'),
                Patch(facecolor='red', alpha=0.2, label='Mode 3')]
ax.legend(handles=legend_modes, loc='upper right', title='Modes')

# Ajout de petites annotations sur le graphique pour identifier les courbes
ax.text(times[-1]*0.7, np.max(Q_abs_full)*0.9, "Q_absorbed", color='black', fontsize=10)
ax.text(times[-1]*0.7, np.max(Q_rel_full)*0.8, "Q_released", color='magenta', fontsize=10)
ax.text(times[-1]*0.7, np.max(Q_loss_full)*0.7, "Q_loss", color='orange', fontsize=10)

plt.tight_layout()
plt.savefig('Q_all_over_time_with_modes_colored.png')
plt.show()