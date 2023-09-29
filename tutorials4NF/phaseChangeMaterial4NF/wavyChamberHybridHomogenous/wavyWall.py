#------------------------------- nanoFluid4Foam project -------------------------------#
#Author
    #Ehsan Golab, SUT. All rights reserved.
    #Ehsan1996Golab@gmail.com

#--------------------------------------------------------------------------------------#
import  math
from numpy import *
import matplotlib.pyplot as plt
import array as arr

length = 1.0
waveNum = 3
amplitude = 0.1
pointNum = 81
deltaX = length/(pointNum-1)
x = arr.array('d',[])
y = arr.array('d',[])
z = 0.1

for i in range(pointNum):
    x.append(i*deltaX)
    y.append(amplitude*math.sin((2*waveNum*math.pi*x[i])/length))


with open("points0_3.h", "w") as points:
    for i in range(len(x)):
	points.write("(\t")
	points.write(str(y[i]))
	points.write("\t")
	points.write(str(x[i]))
	points.write("\t")
	points.write(str(-z))
	points.write("\t)")
	points.write("\n")

with open("points4_7.h", "w") as points:
    for i in range(len(x)):
	points.write("(\t")
	points.write(str(y[i]))
	points.write("\t")
	points.write(str(x[i]))
	points.write("\t")
	points.write(str(z))
	points.write("\t)")
	points.write("\n")

plt.plot(y, x)
plt.xlabel('x')
plt.ylabel('y')
#plt.show()
