#!/bin/bash

NITER=10000
parallel -j1 "./fw.x {1} $NITER > fw_{1}.log && mv phi.dat fw_{1}.dat" ::: {-16..0}
