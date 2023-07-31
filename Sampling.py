#!/usr/bin/env python3

import sys, os, pandas as pnd, pyDOE, numpy as np

NSAM=1000
SMAX=4
NVAR = 5
DATA = np.zeros((NSAM,7))

# Perform Latin hypercube sampling (LHS) or simple Monte Carlo (MC) sampling
do_lhs = True

if do_lhs :
    PARAM = (2 * pyDOE.lhs(NVAR-1, NSAM) - 1.0)*0.03*SMAX
    D4    = (2 * pyDOE.lhs(1, NSAM) - 1.0)*0.01*SMAX
else :
    PARAM = np.random.normal(size=(NSAM, NVAR-1))*0.03
    D4 = np.random.normal(size=(NSAM, 1))*0.01

PARAMS=np.append(PARAM,D4,axis=1)

# Create a DataFrame from the simulation results
DF = pnd.DataFrame(PARAMS, columns=("D0", "D1", "D2", "D3", "D4"))
DF.to_pickle("sample.out")