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
def polyfunc(T, a0, a1, a2, a3):
    return a0 + a1*T + a2*T**2 + a3*T**3

#Data for thermophysical properties (rho[kg/m^3], Cp[j/kg.K], kappa[W/m.K]) of therminol XP for 244<T<600 K
wb = xlrd.open_workbook('materialsPropT.xlsx')
sh = wb.sheet_by_index(1)
nrows = sh.nrows
T = np.zeros(nrows-1)
rho = np.zeros(nrows-1)
Cp = np.zeros(nrows-1)
kappa = np.zeros(nrows-1)

for r in range (nrows-1):
    T[r] = sh.cell_value(r+1,0)
    rho[r] = sh.cell_value(r+1,1)
    Cp[r] = sh.cell_value(r+1,2)
    kappa[r] = sh.cell_value(r+1,3)

#Curve fitting
f = kappa
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
