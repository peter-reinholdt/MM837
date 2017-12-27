#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

base = "output.dat.conf"


plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
data = np.loadtxt(base+"0", dtype=np.float64)
x = np.cos(data)
y = np.sin(data)
m = ax.matshow(data, vmin=0.0, vmax=2*np.pi, cmap=plt.cm.hsv)
q = ax.quiver(x,y)
#plt.colorbar()

i=0
while True:
    data = np.loadtxt(base+str(i), dtype=np.float64)
    x = np.cos(data)
    y = np.sin(data)
    m.set_data(data)
    q.U = x.flat[:]
    q.V = y.flat[:]
    plt.pause(0.001)
    i += 100
