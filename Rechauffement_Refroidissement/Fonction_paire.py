import numpy as np
import matplotlib.pyplot as plt 


# On commence par lire la donnée

box = np.array([10,10])
input = "res/EchauffementL10.xyz"
output = "res/paire_EchauffementL10_10pt.log"
nk = 500
nstep = 10

pos = []

with open(input,"r") as file :
    data = file.readlines()

i = -2
for lines in data :
    if i == -2 :
        n = int(lines.strip())
        pos_temp = np.zeros((n,2))
    if i >= 0 and len(lines.split()) == 3:
        line = lines.split() 
        elmt = line[0]
        pos_temp[i] = [float(line[1]),float(line[2])]
    if i >= 0 and len(lines.split()) != 3:
        i = -2
        pos.append(np.copy(pos_temp))
        # print((pos_temp))
    i += 1

print(f"n = {n}")


print(len(pos))


dr = box[0]/(2*nk)
r = np.linspace(0,nk*dr,nk)
h = np.zeros(nk)
g = np.zeros(nk)
rho = n/np.prod(box)
const = np.pi*rho


# =================================================================================================================================================================
#       Chaud
# =================================================================================================================================================================

print("Calculating for high temperatures ...")


h0 = np.zeros(nk)
g0 = np.zeros(nk)
# Calcul du h
for m in range(1090,1100) :
# for m in range(90,100) :
    for i in range(0,n) :
        for j in range(i+1,n) :
            rij = pos[m][i] - pos[m][j]
            rij -= np.round(rij/box)*box #Permet de respecter les pbc
            rij_sq = 0
            for elmt in rij :
                rij_sq += elmt**2
            k = int(np.floor(np.sqrt(rij_sq)/dr) + 1) #On regarde dans quelle intervalle de distance la particule est
            if (k < nk) :
                h0[k] += 2

#normalisation pour g
for k in range(0,nk) :
    g0[k] = h0[k]/(n*nstep)
    r_low = (k-1)*dr
    r_high = r_low+dr 
    h_id = 2*np.pi*(r_high**2 - r_low**2)*rho #Aire de la surface entre r et r+dr
    g0[k] = g0[k]/h_id

for k in range(nk):
    g0[k] = h0[k] / (n*nstep)
    r_low  = k * dr
    r_high = r_low + dr
    h_id = const * (r_high**2 - r_low**2)
    g0[k] = g0[k] / h_id


plt.plot(r,g0,label = "kT = 10")

# =================================================================================================================================================================
#       Froid
# =================================================================================================================================================================

print("Calculating for low temperatures ...")


h0 = np.zeros(nk)
g0 = np.zeros(nk)
# Calcul du h
for m in range(90,100) :
# for m in range(1090,1100) :
    for i in range(0,n) :
        for j in range(i+1,n) :
            rij = pos[m][i] - pos[m][j]
            rij -= np.round(rij/box)*box #Permet de respecter les pbc
            rij_sq = 0
            for elmt in rij :
                rij_sq += elmt**2
            k = int(np.floor(np.sqrt(rij_sq)/dr) + 1) #On regarde dans quelle intervalle de distance la particule est
            if (k < nk) :
                h0[k] += 2

#normalisation pour g
for k in range(0,nk) :
    g0[k] = h0[k]/(n*nstep)
    r_low = (k-1)*dr
    r_high = r_low+dr 
    h_id = const * (r_high**2 - r_low**2) #Aire de la surface entre r et r+dr
    g0[k] = g0[k]/h_id

plt.plot(r,g0,label = "kT = 0.01")
plt.xlabel("r (\u03c3)")
plt.ylabel("g(r)")
plt.title("Fonction de distribution de paire pour L = 10\u03c3")
plt.legend()
plt.show()