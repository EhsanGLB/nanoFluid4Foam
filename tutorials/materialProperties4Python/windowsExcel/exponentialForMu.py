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
import xlrd


#Functions
def expfunc(T, a, b, c):
    return a * np.exp(-b * T) + c

def EXPFUNC(TT, A, B, C):
    return A * np.exp(-B * TT) + C

def expfuncsym(T, a, b, c):
    return a * sym.exp(-b * T) + c

#Data for thermophysical properties (rho[kg/m^3], Cp[j/kg.K], kappa[W/m.K]) of therminol XP for 244<T<600 K
wb = xlrd.open_workbook('materialsPropT.xlsx')
sh = wb.sheet_by_index(1)
nrows = sh.nrows
T = np.zeros(nrows-1)
mu = np.zeros(nrows-1)

for r in range (nrows-1):
    T[r] = sh.cell_value(r+1,0)
    mu[r] = sh.cell_value(r+1,4)
    
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

