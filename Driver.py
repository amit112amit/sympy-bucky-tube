#!/usr/bin/env python

import sympy as sp
import BuckyTube as b

# Make flat plate points
W = 9
H = 3

# Symbolic computations
platePoints = b.generateRollableFlatPlate( W, H )
cylPoints = b.generateCylinderFromPlate( platePoints )

# Calculate the normals for all the points
R = (W-1)/(2*sp.pi)
normals = []
for pt in cylPoints:
    x,_,z = pt
    normals.append( [ sp.simplify( x/R ), 0, sp.simplify( z/R ) ] )

# Get bonds connectivity for energy calculations
cylinder = b.generateCylinderVTK( W, H, platePoints, cylPoints,
                                writeFile=True)
connectivity = b.makeListFromCellArray( cylinder )

# Iterate over all the bonds and calculate the circularity and normality energy
circEn = 0
normEn = 0
for bond in connectivity:
    i,j = bond
    pti = cylPoints[i]
    ptj = cylPoints[j]
    ni = normals[i]
    nj = normals[j]
    m =  sum( [ (nj[z] - ni[z])**2 for z in range(3) ] )
    r =  sum( [ (ptj[z] - pti[z])**2 for z in range(3) ] )
    circTerm = sum( [ (ni[z] + nj[z])*(ptj[z]-pti[z]) for z in range(3) ] )
    kernel = sp.exp( -r*sp.Rational(1,2) )
    circEn += kernel*circTerm**2/r
    normEn += kernel*m

K = sp.Symbol( 'kappa' )
circEn = K*circEn
normEn = K*sp.simplify( normEn )
