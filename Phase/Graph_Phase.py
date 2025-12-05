import numpy as np 
import matplotlib.pyplot as plt
rho = []
V = 8.77*8.77
E = []
with open("res/out_Phase_N.log","r") as file :
    i = 0
    for lines in file :
        if i < 3 :
            i += 1
        else :
            Line = lines.split()
            rho.append(float(Line[3])/float(Line[4]))
            # print(rho)
            E.append(float(Line[2]))

plt.plot(rho,E,label="N Variable")
rho = []
V = 8.77*8.77
E = []

with open("res/out_Phase_V.log","r") as file :
    i = 0
    for lines in file :
        if i < 3 :
            i += 1
        else :
            Line = lines.split()
            rho.append(float(Line[3])/float(Line[4]))
            # print(rho)
            E.append(float(Line[2]))

plt.plot(rho,E,label="V Variable")
plt.legend()
plt.xlabel("\u03c1 (\u03c3"+"$^{-2}$)")
plt.ylabel("Energie (\u03B5)")

plt.show()
