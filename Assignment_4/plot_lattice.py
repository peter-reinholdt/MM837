#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

base = "output.dat.conf"


plt.ion()
data = np.loadtxt(base+"0", dtype=np.float64)
m = plt.matshow(data, vmin=0.0, vmax=2*np.pi, cmap=plt.cm.jet)
plt.colorbar()

i=0
while True:
    m.set_data(np.loadtxt(base+str(i), dtype=np.float64))
    plt.pause(0.001)
    i += 100
