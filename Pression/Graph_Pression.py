import numpy as np 
import matplotlib.pyplot as plt


kT = 1
N = 100
length = 101

V = np.zeros(length)
P = np.zeros(length)
P_GP = np.zeros(length)

with open("Pressure.log","r") as file :
    i = -1
    for lines in file :
        if i >= 0 :
            Line = lines.split()
            V[i] = Line[0]
            P[i] = Line[1]
            P_GP[i] = N*kT/V[i]
        i += 1


plt.plot(V,P,label="kT = 1")
plt.plot(V,P_GP,label="GP")
plt.legend()
plt.title("Pression selon le Volume")
plt.xlabel("Volume (\u03c3²)")
plt.ylabel("Pression (\u03b5/\u03c3²)")
plt.show()
