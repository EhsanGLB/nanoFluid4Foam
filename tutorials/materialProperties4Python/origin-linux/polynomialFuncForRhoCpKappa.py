#------------------------------- nanoFluid4Foam project -------------------------------#
#Author
    #Ehsan Golab, SUT. All rights reserved.
    #Ehsan1996Golab@gmail.com

#--------------------------------------------------------------------------------------#
import numpy as np
import sympy as sym
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from decimal import Decimal


#Functions
def polyfunc(T, a0, a1, a2, a3):
    return a0 + a1*T + a2*T**2 + a3*T**3


#Data for thermophysical properties (rho[kg/m^3], Cp[j/kg.K], kappa[W/m.K]) of therminol XP for 244<T<600 K
T = np.array([244, 255, 266, 277, 289, 300, 311, 322, 333, 344, 355, 366, 377, 389, 400, 411, 422, 433, 444, 455, 466, 477, 489, 500, 511, 522 ,533, 544, 555 ,566, 577, 589, 600])
rho = np.array([909, 902, 895, 888, 881, 874, 867, 860, 853, 846, 839, 832, 824, 817, 810, 803, 796, 788, 781, 773, 766, 758, 751, 743, 735, 728, 720, 712, 703, 695, 687, 678, 669])
Cp = np.array([1570, 1630, 1680, 1740, 1790, 1850, 1900, 1950, 2000, 2060, 2110, 2160, 2200, 2250, 2300, 2350, 2390, 2440, 2480, 2530, 2570, 2610, 2660, 2700, 2740, 2780, 2820, 2860, 2890, 2930, 2970, 3000, 3040])
kappa = np.array([0.1255, 0.1252, 0.1249, 0.1245, 0.1241, 0.1237, 0.1232, 0.1227, 0.1222, 0.1216, 0.121, 0.1204, 0.1198, 0.1191, 0.1184, 0.1176, 0.1169, 0.1161, 0.1153, 0.1144, 0.1135, 0.1126, 0.1116, 0.1107, 0.1097, 0.1086, 0.1076, 0.1065, 0.1053, 0.1042, 0.103, 0.1018, 0.1005])

#Curve fitting
f = rho
popt, pcov = curve_fit(polyfunc, T, f)

#To show function in above of plot
Ts = sym.Symbol('T')
tex = sym.latex(polyfunc(Ts,*popt)).replace('$', '')
print("f(T) = a0 + a1*T + a2*T^2 + a3*T^3")
print("a0: ", '%.4E' % Decimal(popt[0]))
print("a1: ", '%.4E' % Decimal(popt[1]))
print("a2: ", '%.4E' % Decimal(popt[2]))
print("a3: ", '%.4E' % Decimal(popt[3]))


#To show 
plt.plot(T, f, 'ko', label="Original Data")
plt.plot(T, polyfunc(T, *popt), label="Fitted Curve")
plt.title(r'$f(T)= %s$' %(tex),fontsize=8)
plt.legend()
plt.show()
