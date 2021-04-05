import pip
import numpy as np
from mip import Model, xsum, maximize, BINARY
pip.main(['install', 'mip'])

nt = 18
I  = range(nt) # from 0 to 17
k  = 2 # maximize for team k, we start at 0!!

m = Model("Sports League")

xw = []
xd = []
xl = []

for i in I:
    xw.append([])
    xd.append([])
    xl.append([])
    for j in I:
        if i!=j:
            xw[i].append(m.add_var(var_type=BINARY))
            xd[i].append(m.add_var(var_type=BINARY))
            xl[i].append(m.add_var(var_type=BINARY))
            m += (xw[i][j] + xd[i][j] + xl[i][j]) == 1
        else:
            xw[i].append([])
            xd[i].append([])
            xl[i].append([])

for t in range(nt-1):
    m +=    xsum( 3*xl[i][t]   + xd[i][t]   +
                  3*xw[t][i]   + xd[t][i]  for i in I if i!=t) <=\
            xsum( 3*xl[i][t+1] + xd[i][t+1] +
                  3*xw[t+1][i] + xd[t+1][i]  for i in I if i!=(t+1))

m.objective = maximize(xsum( 3*xl[i][k] + xd[i][k] +
                             3*xw[k][i] + xd[k][i]  for i in I if i!=k))
m.optimize()

pk = np.sum([3*int(xl[i][k].x) + int(xd[i][k].x) +
             3*int(xw[k][i].x) + int(xd[k][i].x) for i in I if i!=k])

print("Team k={} has {} points and will still relegate.".format(k,pk)+\
      "Therefore a team is safe if it has at least {} points".format(pk+1))