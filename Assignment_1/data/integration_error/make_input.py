#!/usr/bin/env python

import numpy as np

totaltime = 10
timesteps = np.logspace(-4,0,200)

for i, ts in enumerate(timesteps):
    with open("{}.inp".format(i), "w") as f:
        out = ""
        out += "[[integrators]]\n"
        out += "dt = {}\n".format(ts) 
        out += "nsteps = {}\n".format(int(totaltime/ts)) 
        out += 'integrator = "velocityVerlet"\n'
        out += "[[system]]\n"
        out += "Nparticles = 21\n"
        out += "delta = 0.33\n"
        out += "k1 = 1.0\n"
        out += "k2 = 1.0\n"
        out += "k3 = 1.0\n"
        out += "[[output]]\n"
        out += "ifreqout = 1\n"
        out += 'outfile  = "{}.dat"\n'.format(i)
        f.write(out)
