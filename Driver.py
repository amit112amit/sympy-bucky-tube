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
#circEn = 0
normEn = 0
for bond in connectivity:
    i,j = bond
    pti = cylPoints[i]
    ptj = cylPoints[j]
    ni = normals[i]
    nj = normals[j]
    m =  (nj[0] - ni[0])**2 + (nj[1] - ni[1])**2 + (nj[2] - ni[2])**2
    #p = [ nj[z] + ni[z] for z in range(3) ]
    rij = [ ptj[z] - pti[z] for z in range(3) ]
    r = (pti[0] - ptj[0])**2 + (pti[1] - ptj[1])**2 + (pti[2] - ptj[2])**2
    kernel = sp.exp( -r*sp.Rational(1,2) )
    #circTerm = ( p[0]*rij[0] + p[1]*rij[1] + p[2]*rij[2] )**2/r
    #circEn += circEn + kernel*circTerm
    normEn += normEn + kernel*m

K = sp.Symbol( 'kappa', real=True, positive=True )
#circEn = K*sp.simplify( circEn )
normEn = K*sp.simplify( normEn )
