import numpy as np
import matplotlib.pyplot as plt 

box = np.array([10, 10])
input = "res/RefroidissementL10.xyz"
nk = 500

# Lecture
pos = []
with open(input,"r") as f:
    lines = f.readlines()

i = -2
for L in lines:
    if i == -2:
        n = int(L.strip())
        pos_temp = np.zeros((n, 2))
    if i >= 0 and len(L.split()) == 3:
        x,y = map(float, L.split()[1:])
        pos_temp[i] = [x, y]
    if i >= 0 and len(L.split()) != 3:
        pos.append(pos_temp.copy())
        i = -2
    i += 1

# Paramètres
dr = box[0] / (2 * nk)
r = np.linspace(0, nk*dr, nk)
rho = n / np.prod(box)

# Histogramme 2D
h = np.zeros(nk)
frames = range(90, 100)   # 10 frames

for m in frames:
    for i in range(n):
        for j in range(i+1, n):
            rij = pos[m][i] - pos[m][j]
            rij -= np.round(rij / box) * box
            d = np.linalg.norm(rij)
            if d < box[0]/2:
                k = int(d / dr)
                if k < nk:
                    h[k] += 2

# Normalisation 2D
g = np.zeros(nk)
for k in range(nk):
    r_mid = (k + 0.5) * dr
    shell = 2 * np.pi * r_mid * dr * rho
    g[k] = (h[k] / len(frames)) / (n * shell)

plt.plot(r, g)
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("RDF en 2D")
plt.show()
