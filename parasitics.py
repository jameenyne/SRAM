#!/usr/bin/env python3

import getopt, numpy as np, pylab as plt, pandas as pnd, os, sys, scipy, pyDOE

def get_nested_cap(W, S, H, T, L) :
    """get the capacitance of a wire of width w spaced by S from wires on the left and right of width W
       with a thickness of T, and a height off the ground plane of H. L is the length. returns fF.
       assumes all dimensions are in microns"""
    C1 = 1.15*(W/H) + 2.8*(T/H)**0.222
    C3 = 2*(0.03*(W/H) + 0.83*(T/H) - 0.07*(T/H)**0.222)*(S/H)**-1.34
    EOX = 3.9 * 8.855e-14
    return EOX*(C1+C3)*(L*1e-6)*1e15

def get_resistance(W, L) :
    """get the resistance of a wire of length L and width W in microns. returns ohms."""
    return 0.05*L/W

print(f"C per segment is {get_nested_cap(0.5,0.5,3.288,0.85,256):.4f}fF {get_resistance(0.5, 2):.4f}")
