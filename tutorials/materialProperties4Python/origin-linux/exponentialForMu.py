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
def expfunc(T, a, b, c):
    return a * np.exp(-b * T) + c

def EXPFUNC(TT, A, B, C):
    return A * np.exp(-B * TT) + C

def expfuncsym(T, a, b, c):
    return a * sym.exp(-b * T) + c

#Data for thermophysical properties (mu[Pa.s]) of therminol XP for 244<T<600 K
T = np.array([244, 255, 266, 277, 289, 300, 311, 322, 333, 344, 355, 366, 377, 389, 400, 411, 422, 433, 444, 455, 466, 477, 489, 500, 511, 522 ,533, 544, 555 ,566, 577, 589, 600])
mu = np.array([6.07, 1.406, 0.433, 0.1661, 0.0755, 0.0392, 0.0227, 0.01423, 0.00956, 0.00679, 0.00504, 0.00388, 0.00307, 0.0025, 0.00207, 0.00175, 0.0015, 0.001301, 0.00114, 0.001008, 0.000898, 0.000805, 0.000726, 0.000658, 0.000599, 0.000547, 0.000502, 0.000461, 0.000425, 0.000393, 0.000363, 0.000337, 0.000313])

#Normalizing data
TT = (T - min(T)) / (max(T) - min(T))
MU = (mu - min(mu)) / (max(mu) - min(mu))

#Curve fitting for normal function
POPT, PCOV = curve_fit(EXPFUNC, TT, MU)
A = POPT[0]
B = POPT[1]
C = POPT[2]

#Convert normal to dimension curve fitting
a = (max(mu) - min(mu)) * A * np.exp(B*min(T)/(max(T) - min(T)))
b = B/(max(T) - min(T))
c = (max(mu) - min(mu)) * C + min(mu)
popt = np.array([a, b, c])

#To show function in above of plot
Ts = sym.Symbol('T')
tex = sym.latex(expfuncsym(Ts,*popt)).replace('$', '')
print("f(T) = a*e^(-b*T) + c")
print("a: ", '%.4E' % Decimal(popt[0]))
print("b: ", '%.4E' % Decimal(popt[1]))
print("c: ", '%.4E' % Decimal(popt[2]))

#To show 
plt.plot(T, mu, 'ko', label="Original Data")
plt.plot(T, expfunc(T, *popt), label="Fitted Curve")
plt.title(r'$f(T)= %s$' %(tex),fontsize=8)
plt.legend()
plt.show()

