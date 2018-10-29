#Mesh Pseudo Curvature Analysis v0.31 - 17/12/2013
#by Mirco Becker - mirco.becker@informance-design.com

import math
import random
import rhinoscriptsyntax as rs
import Rhino
import Rhino.Geometry as rg

obj = Mesh

commonMesh = rg.Mesh()
commonMesh = obj
minRad = []
maxRad = []
mean = []
gausian = []
minDir = []
maxDir = []
geo = []


vertices = commonMesh.Vertices
normals = Normals

for i in range(vertices.Count):

    neigbours = vertices.GetConnectedVertices(i)
    conCount = neigbours.Count

    minC = 10000.0
    maxC = -10000.0

    for j in range(neigbours.Count):

        vVec3d = rg.Vector3d(rg.Point3f.Subtract(vertices[i],vertices[neigbours[j]]))
        edgeLength = vVec3d.Length
        vVec3d.Unitize
        vNor3d = rg.Vector3d(normals[i])
        vAngle = rg.Vector3d.VectorAngle(vVec3d,vNor3d)
        vAngle = vAngle - math.pi * 0.5

        if edgeLength > 0.0:

            if not Rhino.RhinoMath.IsValidDouble(vAngle):
                vAngle  = 0.0
                curvature = 0.0
                
            if vAngle < x and vAngle > (-1.0) * x:
                curvature = 0.0
            else:
                curveRad = edgeLength / math.sin(vAngle)
                curvInv = 1/curveRad

            if curvInv< -1.0:
                curvInv = -1.0

            if curvInv < minC:
                minC = curvInv
                minD = vVec3d
            if curvInv > maxC:
                maxC = curvInv
                maxD = vVec3d

    minRad.append(maxC)
    maxRad.append(minC)
    mean.append((minC+maxC)*0.5)
    gausian.append(minC*maxC)
    minDir.append(minD)
    maxDir.append(maxD)
