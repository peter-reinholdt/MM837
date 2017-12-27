#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]

data = np.loadtxt(filename)

plt.plot(data)
plt.show()
