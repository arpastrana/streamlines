__name__ = "Fungi 6.0"
__author__ = "Rafael Pastrana"
__version__ = "2018.05.07"

import rhinoscriptsyntax as rs


class Fungi(object):

    def __init__(self, curves, minSpacing=200., maxSpacing=400., autoGrade=True):

        self.initialCurves = [curves[0], curves[1]]
        self.seam = rs.AddTweenCurves(curves[0], curves[1])
        self.minSpacing = minSpacing
        self.maxSpacing = maxSpacing
        self.autoGrade = autoGrade

        self.points = []
        self.curves = []
        self.intermediateLines = []
        self.tweenedCurves = []
        self.outputCurves = []

        self.veinsPoints = []
        self.veins = []
        self.thresholds = []
        self.gradingParameters = []

        self.veinsPointsDict = {}
        self.veinsDict = {}
        self.thresholdsDict = {}
        self.clustersDict = {}
        self.connectionsDict = {}


    def recursiveVeinGrowth(self, divisions, iterations, parameters=[]):
        count = 0
        self.gradedVeins(divisions, parameters)
        print(self.thresholdsDict)
        #print(self.connectionsDict)

        while count < iterations:
            # iterate over side A and side B
            for key, thresholds in self.thresholdsDict.items():
                for index, threshold in thresholds.items():
                    #print(threshold)
                    try:
                        curves = self.clustersDict.get(key).get(index)
                    except:
                        break

                    # additional 10 for geometrical stability problems
                    newCurves = self.tweenConnectorCurve(curves, threshold + 0.01)
                    newCurves = self.mergeCurves(curves, newCurves)

                    # replace entry in initial dictionary list
                    try:
                        self.clustersDict[key][index] = newCurves
                        self.outputCurves.extend(newCurves)
                    except:
                        break

                    # increase the count
                    count += 1


    def recursiveCustomVeinGrowth(self, curves, iterations, threshold=320.):
        count = 0
        self.makeCustomVeins(curves, threshold)


        while count < iterations:
            for key, thresholds in self.thresholdsDict.items():
                for index, threshold in thresholds.items():
                    try:
                        curves = self.clustersDict.get(key).get(index)
                    except:
                        break


                    newCurves = self.tweenConnectorCurve(curves, threshold + 0.0)
                    newCurves = self.mergeCurves(curves, newCurves)

                    # replace entry in initial dictionary list
                    try:
                        self.clustersDict[key][index] = newCurves
                        #self.outputCurves.extend(newCurves)
                    except:
                        break

                    # increase the count
                    count += 1


    def retrieveCurvesfromClusters(self):
        outputCurves = []
        #print(self.clustersDict)
        for i in range(len(self.clustersDict)):
            print("i value is: {}".format(i))
            for j in range(len(self.clustersDict[str(i)])):
                print("j value is: {}".format(j))
                outputCluster = self.clustersDict.get(str(i)).get(j)
                print(outputCluster)
                outputCurves.append(outputCluster)


        # for key, cluster in self.clustersDict.items():
        #     for index, curves in cluster.items():



        self.outputCurves = outputCurves


    def makeCustomVeins(self, curves, threshold):
        self.veins = curves
        self.thresholds = [threshold for x in range(len(self.veins)-1)] # -1?
        for idx, curve in enumerate(self.veins):
            if idx < (len(self.veins) - 1):
                cluster = [self.veins[idx], self.veins[idx+1]]
                self.thresholdsDict[str(idx)] = {k:v for k,v in enumerate([self.thresholds[idx]])}
                self.clustersDict[str(idx)] = {k:v for k,v in enumerate([cluster])}


    def createThresholdClusters(self, idx, veins):
        if len(self.thresholds) >= 1:
            clustersDict = {}
            for j, threshold in enumerate(self.thresholds):
                if j < len(self.thresholds):
                    clustersDict[j] = [veins[j], veins[j+1]]
            self.clustersDict[str(idx)] = clustersDict


    def findCustomRoot(self, curve, seam):
        roots = []
        boundaryPoints = self.getCurveEndPoints(curve)
        seamPoints = self.getCurveEndPoints(seam)
        rootPoints = zip(seamPoints, boundaryPoints)

        for rootPointsTuple in rootPoints:
            root = rs.AddLine(rootPointsTuple[0], rootPointsTuple[1])
            roots.append(root)

        roots = self.sortCurvesByLength(roots)
        return roots[0]


    def tweenConnectorCurve(self, curves, threshold):
        """ Iterative finding intermediate points between 2 curves.
        A region for comparison is created according to the smallest curve.
        If distance is lower than threshold, create a point.
        """

        tweenCurves = []
        for index,curve in enumerate(curves):
            if index < (len(curves) - 1):

                # extract subCurve according to smaller one
                sortedCurves = self.sortCurvesByLength([curves[index],
                                                       curves[index+1]
                                                       ]
                                                       )
                # sort curves
                smallCurve = sortedCurves[0]
                largeCurve = sortedCurves[1]
                # subCurve = largeCurve

                # evaluation parameters for shattering
                shatterPoints = self.getCurveEndPoints(smallCurve)
                shatterParam = self.getEvaluationParameters(largeCurve,
                                                            shatterPoints
                                                            )
                #print('my shatter param are {}'.format(shatterParam))
                #print('my shatter points are {}'.format(shatterPoints))

                #find closes point among control points
                newShatterPoints = []
                for shatterPoint in shatterPoints:
                    cloudIndex = rs.PointArrayClosestPoint(rs.CurvePoints(largeCurve), shatterPoint)
                    newShatterPoint = rs.CurvePoints(largeCurve)[cloudIndex]
                    newShatterPoints.append(newShatterPoint)

                shatterParam = self.getEvaluationParameters(largeCurve,
                                            newShatterPoints
                                            )



                #shatterParam = list(map(lambda x: rs.PointClosestObject(x, rs.CurvePoints(largeCurve)),
                #                        shatterPoints))

                #shatterParam = list(map(lambda x: x[1], shatterParam))
                #rs.PointClosestObject(point, object_ids)
                #PointArrayClosestPoint(points, test_point)

                try:
                    subCurve = rs.AddSubCrv(largeCurve, shatterParam[0], shatterParam[1])
                except:
                    subCurve = largeCurve

                # extract control points
                thisPoints = rs.CurvePoints(smallCurve)
                nextPoints = rs.CurvePoints(subCurve)
                pointTuples = zip(thisPoints, nextPoints)

                # find distances between points' pairs
                newPoints = []
                newLines = []

                for idx, pointTuple in enumerate(pointTuples):
                    distance = rs.Distance(pointTuple[0], pointTuple[1])

                    if distance >  threshold:
                        # create average point
                        newLine = rs.AddLine(pointTuple[0], pointTuple[1])
                        parameter = self.getCurveParameter(newLine, 0.5)
                        newPoint = rs.EvaluateCurve(newLine,parameter)
                        newPoints.append(newPoint)


                        self.points.append(newPoint)
                        newLine = rs.AddPolyline([pointTuple[0], newPoint, pointTuple[1]])
                        newLines.append((idx, newLine))
                        self.curves.append(newLine)


                if len(newPoints) > 1:
                    newCurve = rs.AddPolyline(newPoints)
                    tweenCurves.append(newCurve)

        return tweenCurves


    def tweenCurve(self, curves, divisions, threshold):
        """ Iterative finding intermediate points between 2 curves.
        A region for comparison is created according to the smallest curve.
        If distance is lower than threshold, create a point.
        """

        tweenCurves = []
        for index,curve in enumerate(curves):
            if index < (len(curves) - 1):

                # extract subCurve according to smaller one
                sortedCurves = self.sortCurvesByLength([curves[index],
                                                       curves[index+1]
                                                       ]
                                                       )
                # sort curves
                smallCurve = sortedCurves[0]
                largeCurve = sortedCurves[1]

                # evaluation parameters for shattering
                shatterPoints = self.getCurveEndPoints(smallCurve)
                shatterParam = self.getEvaluationParameters(largeCurve,
                                                            shatterPoints
                                                            )

                try:
                    subCurve = rs.AddSubCrv(largeCurve, shatterParam[0], shatterParam[1])
                except:
                    subCurve = largeCurve

                # make points
                thisPoints = rs.DivideCurve(smallCurve, divisions, False, True)
                nextPoints = rs.DivideCurve(subCurve, divisions, False, True)
                pointTuples = zip(thisPoints, nextPoints)

                # find distances between points' pairs
                newPoints = []
                for pointTuple in pointTuples:
                    distance = rs.Distance(pointTuple[0], pointTuple[1])

                    if distance > threshold:
                        # create average point
                        newLine = rs.AddLine(pointTuple[0], pointTuple[1])
                        adaptiveParameter = self.getCurveParameter(newLine, 0.5 * 1.0)
                        newPoint = rs.EvaluateCurve(newLine,adaptiveParameter)

                        newPoints.append(newPoint)
                        self.points.append(newPoint)



                if len(newPoints) > 1:
                    newCurve = rs.AddInterpCurve(newPoints,
                                                 degree=3,
                                                 knotstyle=0,
                                                 start_tangent=None,
                                                 end_tangent=None
                                                 )

                    tweenCurves.append(newCurve)
        return tweenCurves


    def gradedVeins(self, divisions, parameters):
        for idx, curve in enumerate(self.initialCurves):

            veinsPoints = []
            veins = {}
            growthRoot = self.findRoot(curve)
            intermediateLines = self.createIntermediateLines([self.seam, curve],
                                                             divisions
                                                             )

            if self.autoGrade == True:
                self.gradingParameters, self.thresholds = self.getLinearGradingParameters(growthRoot)
            else:
                self.setNewGradingParameters(parameters)
                self.thresholds = self.getLengthThresholds(growthRoot, self.gradingParameters)

            veinsPoints = self.findVeinsPoints(intermediateLines)
            veins = self.createVeinCurves(veinsPoints, veins)


            # create dictionaries
            self.veinsPointsDict[str(idx)] = veinsPoints # is this really necessary?
            self.veinsDict[str(idx)] = veins # is this really necessary?

            self.intermediateLines.extend(intermediateLines)
            self.thresholdsDict[str(idx)] = {k:v for k,v in enumerate(self.thresholds)}
            self.connectionsDict[str(idx)] = {k:v for k,v in enumerate(intermediateLines)}
            self.createThresholdClusters(idx, veins)


    def mergeCurves(self, curves, newCurves):
        """ It weaves curves and tweened curves lists in a 0-1 manner.
            It updates self.curves and clears self.tweenedCurves.
        """

        mergedCurves = self.weaveLists(curves, newCurves)
        newCurves[:] = []
        return mergedCurves


    def setNewSeam(self,curve):
        self.seam = curve


    def setNewGradingParameters(self, parameters):
        """Replaces with a list of at least 2 parameters.
            Always including 0 and 1.
        """
        self.gradingParameters = parameters


    def findVeinsPoints(self, lines):
        points = []
        for line in lines:
            gradedPoints = []

            for gradingParameter in self.gradingParameters:
                param = self.getCurveParameter(line, gradingParameter)
                gradedPoint = rs.EvaluateCurve(line, param)
                gradedPoints.append(gradedPoint)

            points.append(gradedPoints)
            self.veinsPoints.append(gradedPoints)
        return points


    def findRoot(self, curve):
        roots = []
        boundaryPoints = self.getCurveEndPoints(curve)
        seamPoints = self.getCurveEndPoints(self.seam)
        rootPoints = zip(seamPoints, boundaryPoints)

        for rootPointsTuple in rootPoints:
            root = rs.AddLine(rootPointsTuple[0], rootPointsTuple[1])
            roots.append(root)

        roots = self.sortCurvesByLength(roots)
        return roots[0]


    def createVeinCurves(self, points, veins):
        for index in range(len(self.gradingParameters)):
            tempPoints = []
            for gradedPoints in points:
                tempPoints.append(gradedPoints[index])
            vein = rs.AddPolyline(tempPoints)
            veins[index] = vein
            self.veins.append(vein)
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


    def sortCurvesByLength(self, curvesPair):
        sortedCurves = sorted(curvesPair, key = lambda x: rs.CurveLength(x))
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
