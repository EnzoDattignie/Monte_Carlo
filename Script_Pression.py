import os
import numpy as np 
import matplotlib.pyplot as plt

box_size = np.linspace(10,20,101)
P = np.zeros(len(box_size))
lmax = np.linspace(0.03,0.3,101)

V = []
P = []
P_temp = np.zeros(1001)
lmax = lmax*4
print(lmax)

for i in range(0,len(box_size)) :
    os.system("./config_cryst_2D.out 10 "+str(box_size[i]))
    os.system("./a.out 22 1000 "+str(lmax[i])+" 1")
    V.append(box_size[i]*box_size[i])
    with open("res/out.log","r") as file :
        i = -2
        for lines in file :
            if i >= 0 :
                P_temp[i] = float(lines.split()[3])
                
            i += 1
    # print(P_temp)
    P.append(np.mean(P_temp[-100:]))

os.system("rm res/Pressure.log")
with open("res/Pressure.log","w") as file:
    file.write("V P\n")
    for i in range(0,len(V)) :
        file.write(f"{V[i]} {P[i]}\n")
kT = 1
P_gp = np.zeros(len(V))
for i in range(0, len(P_gp)) :
    P_gp[i] = 100*kT/V[i]


plt.plot(V,P,label="kT = 1")
plt.plot(V,P_gp,label="Pression GP")
plt.legend()
plt.show()