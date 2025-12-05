import numpy as np 
import matplotlib.pyplot as plt

n_cycle = []
E = []

with open("res/out.log","r") as file :
    i = 0
    for lines in file :
        if i < 3 :
            i += 1
        else :
            Line = lines.split()
            n_cycle.append(float(Line[0]))
            E.append(float(Line[2]))

plt.plot(n_cycle,E)
plt.show()
