#!/usr/bin/env python3

import sympy as sp

sp.init_printing()

def phiN(W,H,K):
    t = sp.pi/W
    N = sp.simplify( 2*(H+1)*W*K*sp.exp( sp.Rational(-1,2) )*(1 - sp.cos(2*t))\
        + 4*H*W*K*sp.exp( sp.Rational(-1,2) )*(1 - sp.cos(t)))
    return N


def kappa(W,H,K):
    out = phiN(W,H,K)/H/sp.pi/sp.sqrt(3)
    return out
