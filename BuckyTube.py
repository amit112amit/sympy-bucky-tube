#!/usr/bin/env python3
"""
This module provides function to create a zig-zag Bucky tube
"""

import sympy as sp
import vtk as v
import numpy as np

"""
The following function creates points lying on a flat plate
that can be rolled up into a cylinder. We assume that the plate
is lying in the x-y plane at a constant value of z = -R
r : lattice spacing for the triangular lattice
W : Width in units of r ( should be odd )
H : Height in units of r
R : Radius of the cylinder formed after rolling-up
To get a closed cylinder R == (W - 1)/(2*pi)
"""
def generateRollableFlatPlate( r, W, H, R ):
    points = []
    z = -R

    for i in range(H):
        y = -i*r*sp.sqrt(3)*sp.Rational(1,2)
        if i%2 == 0:
            x = 0
        else:
            x = r*sp.Rational(1,2)

        for j in range(W):
            p = [ x + j*r, y, z ]
            points.append( p )

    return points

"""
The following function calculates coordinates of a cylinder formed by
rolling up points on a flat plate obtained by 'generateRollableFlatPlate'
function.
"""
def generateCylinderFromPlate( points ):
    R = -points[0][2]
    cylinderPoints = []
    for p in points:
        theta = p[0]/R - sp.pi*sp.Rational(1,2)
        c = [ R*sp.cos( theta ), p[1], R*sp.sin( theta ) ]
        cylinderPoints.append(c)
    return cylinderPoints

"""
The following function makes polydata with triangles from a
point cloud lying on a flat plate generated using generateRollableFlatPlate.
"""
def generateMeshPolyData( points ):
    pts = v.vtkPoints()
    for p in points:
        pts.InsertNextPoint( p )
    pd = v.vtkPolyData()
    pd.SetPoints( pts )
    # Create a IdFilter to preserve original point ids
    idf = v.vtkIdFilter()
    idf.SetIdsArrayName( 'origIds' )
    idf.PointIdsOn()
    idf.SetInputData( pd )
    # Create a 2D Delaunay triangulation
    d2d = v.vtkDelaunay2D()
    d2d.SetInputConnection( idf.GetOutputPort() )
    # Create mesh quality cell data array to filter out
    # boundary triangles
    mq = v.vtkMeshQuality()
    mq.SetInputConnection( d2d.GetOutputPort() )
    mq.SetTriangleQualityMeasureToAspectRatio()
    mq.SaveCellQualityOn()
    mq.Update()
    # Generate the polydata with mesh
    plateTriPoly = mq.GetOutput()
    plateTri = plateTriPoly.GetPolys()
    # Map the connectivity to original point ids
    newTriangles = v.vtkCellArray()
    plateTri.InitTraversal()
    triIds = v.vtkIdList()
    origIds = plateTriPoly.GetPointData().GetArray('origIds')
    aspectRatioArray = plateTriPoly.GetCellData().GetArray( 'Quality' )
    triangleIndex = 0
    while plateTri.GetNextCell( triIds ):
        # Keep only the equilateral triangles to get rid of boundary triangles
        aspectRatio = aspectRatioArray.GetTuple1( triangleIndex )
        triangleIndex += 1
        if aspectRatio < 1.5:
            newTriangles.InsertNextCell( 3 )
            for i in range(3):
                newTriangles.InsertCellPoint( int( origIds.GetTuple1(
                    triIds.GetId(i) ) ) )
    pd.SetPolys( newTriangles )
    return pd

"""
The following function uses the above SymPy function
to calculate numeric values for points on a flat plate and
writes it to a vtkPolyData file for visualization in Paraview.
"""
def generatePlateVTK( r, W, H, R, fileName ):
    # Symbolic computations
    rs,Rs = sp.symbols( 'r,R' )
    platePoints = generateRollableFlatPlate( rs, W, H, R )
    # Numerical computations
    plateFunc = sp.lambdify( (rs, Rs), platePoints )
    realPlatePoints = np.array( plateFunc( r, R ) )
    # Generate connectivity
    pd = generateMeshPolyData( realPlatePoints )
    # Write PolyData to file
    w = v.vtkPolyDataWriter()
    w.SetInputData( pd )
    w.SetFileName( fileName )
    w.Write()
    return pd


"""
The following function uses the above two SymPy functions
to calculate numeric values for points on a cylinder and
writes it to a vtkPolyData file for visualization in Paraview.
"""
def generateCylinderVTK( r, W, H, R, fileName ):
    # Symbolic computations
    rs,Rs = sp.symbols( 'r,R' )
    platePoints = generateRollableFlatPlate( rs, W, H, R )
    cylinderPoints = generateCylinderFromPlate( platePoints )
    # Numerical computations
    plateFunc = sp.lambdify( (rs, Rs), platePoints )
    cylFunc = sp.lambdify( (rs, Rs), cylinderPoints )
    realPlatePoints = np.array( plateFunc( r, R ) )
    realCylPoints = np.array( cylFunc( r, R ) )
    # Generate delaunay triangulation of the flat plate
    plateTri = generateMeshPolyData( realPlatePoints )
    cylPoints = v.vtkPoints()
    for p in realCylPoints:
        cylPoints.InsertNextPoint( p )
    cylPoly = v.vtkPolyData()
    cylPoly.SetPoints( cylPoints )
    cylPoly.SetPolys( plateTri.GetPolys() )
    # We have coincident points at the seam that we need to clean
    cd = v.vtkCleanPolyData()
    cd.SetInputData( cylPoly )
    cd.SetTolerance( 1e-3 )
    cd.PointMergingOn()
    cd.Update()
    # Write the cylinder to a VTK file
    w = v.vtkPolyDataWriter()
    w.SetInputData( cd.GetOutput() )
    w.SetFileName( fileName )
    w.Write()
    return

if __name__ == "__main__":
    """
    Test the functions that we wrote above
    """
    W = 13
    H = 7
    r = 1.0
    R = (W-1)/(2*np.pi)
    generatePlateVTK( r, W, H, R, 'plate.vtk' )
    print('Generating cylinder.')
    generateCylinderVTK( r, W, H, R, 'cyl.vtk' )
