# Import our modules that we are using
import numpy as np
import matplotlib.pyplot as plt

# GENERAL INFO
cover = 2
db = 1.27
db_tie = 0.375
h = 24
b = 24
Eu = 0.003
fc = 5000
fy = 60000
Es = 29000
Ag = h * b
phi2 = 0.8
phimax = 0.65
Type = 1

if Type >= 1:
    theta = 0.8
else:
    theta = 0.85

# Matrix info
M = np.array(
[[1, 1, 1, 1],
 [1, 0, 0, 1],
 [1, 0, 0, 1],
 [1, 1, 1, 1]])

# number of columns
ns = np.size(M,0)


# area of steel
ns1 = np.sum(M,0)  #sum bars per column
ns2 = sum(ns1)     # sum all bars
As = (np.pi * db**2) / 4  # area of 1 bars
Ast = ns2 * As



# B1
if fc <= 4000:
    B1 = 0.85
elif fc >= 8000:
    B1 = 0.65
else:
    B1 = (0.85 - (0.00005 * (fc - 4000)))


# COVER PLUS DELTA
d1 = cover + (db/2) + db_tie
delta = (h - d1*2)/(ns - 1)

i = np.arange(1,ns+1,1)
a = np.arange(1,h,.2)

# d
d = d1 + (delta *(i-1))
dt = np.max(d)
# ct
ct = a / B1

fs = np.arange(i.size * a.size).reshape(i.size, a.size)
for ii in range(i.size):
    for aa in range(a.size):
        es = Eu * (d[ii] - ct[aa]) / ct[aa]
        fs[ii, aa] = np.sign(es) * min(np.abs(Es * es), fy)


# Pn Max

Pnmax = theta * phimax * (0.85 * fc * (Ag - Ast) + (Ast * fy))

#phiA = np.arange(a.size)
Pn = np.arange(a.size)
Mn = np.arange(a.size)
for aa in range(a.size):
    et = Eu * ( dt - ct[aa]) / ct[aa]
    if et > 0.0005:
        phiA = 0.9
    elif et < 0.0002:
        phiA = 0.65
    else:
        phiA = (1.45 + 200 * et) / 3

    pnx = 0
    mnx = 0
    for ii in range(i.size):
        pnx += As * ns1[ii] * fs[ii][aa]
        mnx += As * ns1[ii] * fs[ii][aa] * (d[ii]-h/2)
    Pn[aa] = min((phiA * 0.85 * fc * a[aa] * b - pnx), Pnmax)
    Mn[aa] = phiA * (0.85 * fc * a[aa] * b * (h/2 - a[aa]/2) + mnx)


plt.plot(Pn, Mn)

# naming the x axis
plt.xlabel('Moment (ft*Kips)')
# naming the y axis
plt.ylabel('Axial (kips)')

# giving a title to my graph
plt.title('Column Interaction Diagram')

# Show the plot
plt.show()