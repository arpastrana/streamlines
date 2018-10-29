__name__ = "Tweener 1.0"
__author__ = "Rafael Pastrana"
__version__ = "2018.05.01"

import rhinoscriptsyntax as rs
import math
from ghpythonlib import components as ghcomp



class Tweener(object):

    def __init__(self):
        self.curves = []


    def adaptiveTweenCurves(self, curves, divisions, min, max, numMin=1, numMax=1, data=None):
        """ Insert Docstrings here.
        """
        outputCurves = []
        for index,curve in enumerate(curves):
            print("--------------------------------")
            print("Current Tree is {}".format(index))

            if index < (len(curves) - 1):
                tempCurves = []
                cluster = [curves[index], curves[index+1]]
                seam = rs.AddTweenCurves(curves[index], curves[index+1])

                for j, curve in enumerate(cluster):
                    print("*****************************")
                    print("Current Side is: {}".format(j))
                    growthRoot = self.findRoot(curve, seam)
                    lines = self.createIntermediateLines([curve, seam], divisions)


                    parameters = self.getConstrainedGradingParameters(growthRoot,
                                                                      min,
                                                                      max,
                                                                      numMin,
                                                                      numMax)
                    if data != None:
                        print("Linear mapped parameters are: {}".format(parameters))
                        parameters = ghcomp.Interpolatedata(data, parameters)
                        print("GH interpolated parameters are: {}".format(parameters))
                        print(parameters)

                    veinsPoints = self.findVeinsPoints(lines, parameters)
                    veins = self.createVeinCurves(veinsPoints, parameters)

                    if j == 0:
                        veins = veins[:-1]
                    elif j == 1:
                        veins.reverse()

                    tempCurves.extend(veins)
                outputCurves.append(tempCurves)

        self.curves = outputCurves


    def getConstrainedGradingParameters(self, curve, min, max, numMin, numMax):

        #outputParameters = [0]
        totalLength = rs.CurveLength(curve)
        print("Total Length is {}".format(totalLength))

        initialConstrains = [min for i in range(numMin)]
        lastConstrains = [max for i in range(numMax)]

        initialValue = 0
        for initialConstrain in initialConstrains:
            initialValue += initialConstrain


        lastValue = 0
        for lastConstrain in lastConstrains:
            lastValue += lastConstrain
        tempLastValue = lastValue
        #print("tempLastValue is {}".format(tempLastValue))
        lastValue = totalLength - lastValue


        # if totalLength - initialValue - tempLastValue <= min + max:
        #     print("MAX WAS MODIFIED!")
        #     max = max / 2.0
        #     lastConstrains = [max for i in range(numMax)]
        #
        #     lastValue = 0
        #     for lastConstrain in lastConstrains:
        #         lastValue += lastConstrain
        #     lastValue = totalLength - lastValue



        initialValues = [0] + self.massAddition(initialConstrains)[:-1]
        lastValues = list(map(lambda x: totalLength-x,self.massAddition(lastConstrains)))
        lastValues.reverse()
        lastValues = lastValues[1:] + [totalLength]

        # print("Initial Value is: {}".format(initialValue))
        # print("Initial Const are: {}".format(initialConstrains))
        # print("Initial Values are: {}".format(initialValues))
        # print("Last Value is: {}".format(lastValue))
        # print("Last Values are: {}".format(lastValues))

        if totalLength - initialValue - tempLastValue <= min + max:
            subCurve = curve
        else:
            subCurve = rs.AddSubCrv(curve, initialValue, lastValue) # check order of input

        length = rs.CurveLength(subCurve)
        count = 2
        search = True

        print("Sub Length is {}".format(length))


        while search:
            evaluationRange = range(count)
            contenderLength = 0.0

            for value in evaluationRange:
                remappedValue = self.linearInterpolation(0, value, count-1, min, max)
                contenderLength += remappedValue

            if contenderLength >= length:
                search = False
            else:
                count += 1


        count = count - 1
        print("Count is: {}".format(count))

        if count < 2:
            values = [0, 1]

        else:
            values = range(count)
            #print("Range is: {}".format(values))
            values = list(map(lambda x: self.linearInterpolation(
                                                                 values[0],
                                                                 x,
                                                                 values[-1],
                                                                 min,
                                                                 max
                                                                 ),
                                                                values
                                                              )
                                                          )

            #print("Initial interpolated values are: {}".format(values))
            values = self.massAddition(values)
            #print("Mass added values are: {}".format(values))

            values = list(map(lambda x: x/length, values))
            #print("Normalized values are: {}".format(values))
            values.insert(0, 0.0) # ??
            #print("1-zero inserted values are: {}".format(values))
            values = list(map(lambda x: self.linearInterpolation(
                                                                 values[0],
                                                                 x,
                                                                 values[-1],
                                                                 initialValue,
                                                                 lastValue
                                                                 ),
                                    values
                              )
                          )
            #print("Reconstrained values are: {}".format(values))

            values = initialValues + values + lastValues
        #print("Joint values are: {}".format(values))

            distances = self.getConsecutiveDiferences(values, reverse=True)
            print("Differences are: {}".format(distances))

            values = list(map(lambda x: self.linearInterpolation(
                                                                         values[0],
                                                                         x,
                                                                         values[-1],
                                                                         0.0,
                                                                         1.0
                                                                         ),
                                            values
                                      )
                                  )

        ### Insert Grasshopper interpolate function
        #print("Final interpolated values are: {}".format(values))

        return values


    def cubicFunction(x, s, t, u, v):
        a = s*math.pow((1-x),3)
        b = 3*t*x*math.pow((1-x),2)
        c = 3*t*(1-x)*math.pow(x,2)
        d = v*math.pow(x,3)
        return a + b+ c + d


    def mergeCurves(self, curves, newCurves):
            """ It weaves curves and tweened curves lists in a 0-1 manner.
                It updates self.curves and clears self.tweenedCurves.
            """

            mergedCurves = self.weaveLists(curves, newCurves)
            newCurves[:] = []
            return mergedCurves


    def findVeinsPoints(self, lines, parameters):
        points = []
        for line in lines:

            gradedPoints = []
            for parameter in parameters:
                param = self.getCurveParameter(line, parameter)
                gradedPoint = rs.EvaluateCurve(line, param)
                gradedPoints.append(gradedPoint)

            points.append(gradedPoints)

        return points


    def findRoot(self, curve, seam):
        roots = []
        boundaryPoints = self.getCurveEndPoints(curve)
        seamPoints = self.getCurveEndPoints(seam)
        rootPoints = zip(seamPoints, boundaryPoints)

        for rootPointsTuple in rootPoints:
            root = rs.AddLine(rootPointsTuple[0], rootPointsTuple[1])
            roots.append(root)

        roots = self.sortCurvesByLength(roots)
        return roots[0]


    def createVeinCurves(self, points, parameters):
        veins = []

        for index in range(len(parameters)):
            tempPoints = []

            for gradedPoints in points:
                tempPoints.append(gradedPoints[index])

            vein = rs.AddPolyline(tempPoints)
            veins.append(vein)

        return veins


    def createIntermediateLines(self, curves, divisions):
        pointsA = rs.DivideCurve(curves[0], divisions, False, True)
        pointsB = rs.DivideCurve(curves[1], divisions, False, True)
        pointTuples = zip(pointsA, pointsB)

        intermediateLines = []
        for pointTuple in pointTuples:
            newLine = rs.AddLine(pointTuple[0], pointTuple[1])
            intermediateLines.append(newLine)
        return intermediateLines


    def linearInterpolation(self, x1, x2, x3, y1, y3):
        return ((x2 - x1) * (y3 - y1) / (x3 - x1)) + y1


    def getLinearGradingParameters(self, curve):
        length = rs.CurveLength(curve)
        count = 2
        search = True

        while search:
            evaluationRange = range(count)
            contenderLength = 0.0

            for value in evaluationRange:
                remappedValue = self.linearInterpolation(0, value, count-1, self.minSpacing, self.maxSpacing)
                contenderLength += remappedValue

            if contenderLength >= length:
                search = False
            else:
                count += 1

        count = count - 1

        reparametrizedValues = range(count)
        reparametrizedValues = list(map(lambda x: self.linearInterpolation(
                                                                     reparametrizedValues[0],
                                                                     x,
                                                                     reparametrizedValues[-1],
                                                                     self.minSpacing,
                                                                     self.maxSpacing
                                                                     ),
                                        reparametrizedValues
                                  )
                              )

        reparametrizedValues = self.massAddition(reparametrizedValues)
        reparametrizedValues = list(map(lambda x: x/length, reparametrizedValues))
        reparametrizedValues.insert(0, 0.0)
        reparametrizedValues = list(map(lambda x: self.linearInterpolation(
                                                             reparametrizedValues[0],
                                                             x,
                                                             reparametrizedValues[-1],
                                                             1.0,
                                                             0.0
                                                             ),
                                reparametrizedValues
                          )
                      )


        nrp = self.getLengthThresholds(curve, reparametrizedValues)
        return reparametrizedValues, nrp


    def getConsecutiveDiferences(self, values, reverse=False):
        if reverse == False:
             return [s - t for s, t in zip(values, values[1:])]
        else:
             return [t - s for s, t in zip(values, values[1:])]


    def getLengthThresholds(self, curve, parameters):
        nrp = list(map(lambda x: x * rs.CurveLength(curve), parameters))
        nrp = self.getConsecutiveDiferences(nrp)
        return nrp


    def massAddition(self, values):
        addedValues = []
        addedValue = 0.

        for value in values:
            addedValue += value
            addedValues.append(addedValue)
        return addedValues


    def getCurveParameter(self, curve, parameter):
        if parameter >= 0. and parameter <= 1.:
            return parameter * rs.CurveLength(curve)


    def weaveLists(self,list1,list2):
        weavedList = []
        i = 0
        while i < len(list1):
            weavedList.append(list1[i])
            if i in range(len(list2)):
                weavedList.append(list2[i])
            i += 1
        return weavedList


    def sortCurvesByLength(self, curves):
        sortedCurves = sorted(curves, key = lambda x: rs.CurveLength(x))
        return sortedCurves


    def getCurveEndPoints(self, curve):
        startPoint = rs.CurveStartPoint(curve)
        endPoint = rs.CurveEndPoint(curve)
        return [startPoint, endPoint]


    def getEvaluationParameters(self, curve, points):
        parameters = []
        for point in points:
            parameter = rs.CurveClosestPoint(curve, point, segment_index=-1)
            parameters.append(parameter)
        return parameters


    def getClosestPointParameters(self, curve, points):
        parameters = []
        controlPoints = rs.CurvePoints(curve)
        for point in points:
            index = rs.PointArrayClosestPoint(controlPoints, point)
            pointOnCurve = controlPoints[index]
            parameter = rs.CurveClosestPoint(curve, pointOnCurve, segment_index=-1)
            parameters.append(parameter)
        return parameters


    def shatterCurve(self, curve, evaluationParameters):
        curveSegments = rs.SplitCurve(curve, evaluationParameters)
        return curveSegments
