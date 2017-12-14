#!/usr/bin/env python
"""
This module provides functions to create a zig-zag Bucky tube
"""

import sympy as sp
import vtk as v
import numpy as np
import collections as col

"""
The following function creates points lying on a flat plate
that can be rolled up into a cylinder. We assume that the plate
is lying in the x-y plane at a constant value of z = -R
W : Width in units of lattice spacing ( should be odd )
H : Height in units of lattice spacing
"""
def generateRollableFlatPlate( W, H ):
    R = (W-1)/(2*sp.pi)
    points = []
    z = -R

    for i in range(H):
        y = -i*sp.sqrt(3)*sp.Rational(1,2)
        if i%2 == 0:
            x = 0
        else:
            x = sp.Rational(1,2)

        for j in range(W):
            p = [ x + j, y, z ]
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
        if c not in cylinderPoints:
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
def generatePlateVTK( platePoints, fileName ):
    # Numerical computations
    plateFunc = sp.lambdify( (), platePoints )
    realPlatePoints = np.array( plateFunc() )
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
def generateCylinderVTK( W, H, platePoints, cylPoints, writeFile=False,
                        fileName='cylinder.vtk' ):
    # Numerical computations
    plateFunc = sp.lambdify( (), platePoints )
    cylFunc = sp.lambdify( (), cylPoints )
    realPlatePoints = np.array( plateFunc() )
    realCylPoints = np.array( cylFunc() )
    # Generate delaunay triangulation of the flat plate
    plateTriPoly = generateMeshPolyData( realPlatePoints )
    plateTri = plateTriPoly.GetPolys()
    # Create a mapping from plate point Ids to cylinder ids
    # to handle the overlapping points on the seam
    plateToCylIndexMap = np.array([i for i in range(H*W)])
    B = np.array([i for i in range(H)])
    plateToCylIndexMap = plateToCylIndexMap.reshape( (H,W) )
    B = B.reshape( (H,1) )
    plateToCylIndexMap = plateToCylIndexMap - B
    plateToCylIndexMap[:,W-1] = plateToCylIndexMap[:,0]
    plateToCylIndexMap = plateToCylIndexMap.reshape( (H*W,) )
    # Now create the cylinder PolyData
    cylPoints = v.vtkPoints()
    for p in realCylPoints:
        cylPoints.InsertNextPoint( p )
    cylPoly = v.vtkPolyData()
    cylPoly.SetPoints( cylPoints )
    cylTri = v.vtkCellArray()
    triIds = v.vtkIdList()
    plateTri.InitTraversal()
    while plateTri.GetNextCell( triIds ):
        pt0 = triIds.GetId(0)
        pt1 = triIds.GetId(1)
        pt2 = triIds.GetId(2)
        cylTri.InsertNextCell( 3 )
        cylTri.InsertCellPoint( plateToCylIndexMap[ pt0 ] )
        cylTri.InsertCellPoint( plateToCylIndexMap[ pt1 ] )
        cylTri.InsertCellPoint( plateToCylIndexMap[ pt2 ] )
    cylPoly.SetPolys( cylTri )
    if writeFile:
        # Write the cylinder to a VTK file
        w = v.vtkPolyDataWriter()
        #w.SetInputData( cd.GetOutput() )
        w.SetInputData( cylPoly )
        w.SetFileName( fileName )
        w.Write()
    return cylPoly

# This function extracts from a VTK polydata object a list of lists
# representing the edges of the polydata
def makeListFromCellArray( poly ):
    listOut = []
    idf = v.vtkIdFilter()
    idf.SetInputData( poly )
    idf.SetIdsArrayName( 'OrigIds' )
    idf.PointIdsOn()
    edgeExt = v.vtkExtractEdges()
    edgeExt.SetInputConnection( idf.GetOutputPort() )
    edgeExt.Update()
    edges = edgeExt.GetOutput()
    origIds = edges.GetPointData().GetArray('OrigIds')
    ids = v.vtkIdList()
    cells = edges.GetLines()
    cells.InitTraversal()
    while cells.GetNextCell( ids ):
        connectivity = []
        for i in range( ids.GetNumberOfIds() ):
            connectivity.append( int( origIds.GetTuple1( ids.GetId( i ) ) ) )
        listOut.append( connectivity )
    return listOut

if __name__ == "__main__":

    #Test the functions that we wrote above
    W = 13
    H = 7
    # Symbolic computations
    platePoints = generateRollableFlatPlate( W, H )
    cylinderPoints = generateCylinderFromPlate( platePoints )

    generatePlateVTK( platePoints, 'plate.vtk' )
    print('Generating cylinder.')
    generateCylinderVTK( W, H, platePoints, cylPoints, writeFile=True,
                        fileName='cyl.vtk' )

