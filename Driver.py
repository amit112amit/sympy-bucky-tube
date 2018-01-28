#!/usr/bin/env python

import sympy as sp
import sympy.physics.vector as v
import BuckyTube as b

sp.init_printing()
# Make flat plate points
W = 20
H = 2

# Symbolic computations
platePoints = b.generateRollableFlatPlate( W, H )
cylPoints = b.generateCylinderFromPlate( platePoints )

# Calculate the normals for all the points
R = W/(2*sp.pi)
normals = []
for pt in cylPoints:
    x,_,z = pt
    normals.append( [ sp.simplify( x/R ), 0, sp.simplify( z/R ) ] )

# Get bonds connectivity for energy calculations
cylinder = b.generateCylinderVTK( W, H, platePoints, cylPoints,
                                writeFile=True)
connectivity = b.makeConnectivityList( cylinder )

# Iterate over all the bonds and calculate the circularity and normality energy
N = v.ReferenceFrame('N')
circEn = 0
normEn = 0
for bond in connectivity:
    i,j = bond
    pti = cylPoints[i][0]*N.x + cylPoints[i][1]*N.y + cylPoints[i][2]*N.z
    ptj = cylPoints[j][0]*N.x + cylPoints[j][1]*N.y + cylPoints[j][2]*N.z
    ni = normals[i][0]*N.x + normals[i][1]*N.y + normals[i][2]*N.z
    nj = normals[j][0]*N.x + normals[j][1]*N.y + normals[j][2]*N.z
    m = (nj - ni).dot( (nj - ni) ) 
    circTerm = ((ni + nj).dot(ptj - pti))**2
    kernel = sp.exp( -sp.Rational(1,2) )
    circEn += kernel*circTerm
    normEn += kernel*m

K = sp.Symbol( 'K' )
circEn = K*circEn
circEn = K*sp.simplify( circEn )
normEn = K*sp.simplify( normEn )
