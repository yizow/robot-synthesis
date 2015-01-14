from scipy.optimize import minimize as optimize

import database
import featureVector

def findClosest(testTrace, traces, distanceMetric):
    testTrace = testTrace
    with open('testFeature.txt', 'w') as f:
        print 'Finding Closest'
        div = 32.0
        progress = 0.0
        increment = 0.0
        distances = []
        counter = 0
        for trace in traces:
            if counter >= increment:
                print "%d%%" % (progress*100)
                increment += len(traces)/div
                progress += 1/div
            counter += 1
            dist = distanceMetric(trace, testTrace)
            f.write(str(dist))
            f.write('\n')
            distances += [dist]
        distances = [[distances[index], index] for index in range(len(distances))]

    distances.sort()
    print 'Closest image: %d' % distances[0][1]
    return distances


def optimizeParametersNM(testTrace, index):
    p = parameterLookup[index]

    def optimizingFunction(p):
        trace = []
        inCrank = Beam2(p[0])
        rocker = Beam2(p[1])
        outCrank = Beam2(p[2])
        base = Beam2(p[3])
        PoIOffset = p[4]
        PoIDistance = p[5]
        rocker.PoIOffset = PoIOffset
        rocker.PoIDistance = PoIDistance
        for angle in range(1,int(NUMPOINTS)):
            a = angle*pi/2.0/(NUMPOINTS/4)
            position = construct.buildState(inCrank, rocker, outCrank, base, a)
            if not position:
                continue
            trace += [rocker.PoI()]
        return getDistanceMetric(testTrace, trace)

    return optimize(optimizingFunction, p, method = 'Nelder-Mead')


def optimizeParameters(testTrace, index, parameterLookup, buildType, distanceMetric):
    """ Takes the testTrace and index of a certain trace's parameters in parameterLookup and further optimizes those parameters
    based on just the distance between points. 

        Let p be a vector representing the parameters we are optimizing over
        p = [inCrankLength rockerLength outCrankLength baseLength PoIOffset PoIDistance]
        Subject to the constraints:
            baseLength > inCrankLength
            baseLength > outCrankLength
        Grashoff Conditions
            baseLength+rockerLength-inCrankLength-outCrankLength > 0 
            baseLength-rockerLength-inCrankLength+outCrankLength > 0
            baseLength+rockerLength-inCrankLength+outCrankLength > 0

    Returns the new parameters as a list
    """
    p = parameterLookup[index]
    delta = .2

    print 'Further Optimizing'
    print 'Generating finer traces'
    
    paramsGen = (p[0]-delta, p[0]+delta, p[1]-delta, p[1]+delta, p[2]-delta, p[2]+delta, p[3]-delta, p[3]+delta, delta)
    offsetGen = (p[4]-delta/3.0, p[4]+delta/3.0, p[5], p[5]+delta, delta)
    fineMechanisms, finePoI, fineParam = [], [], []

    database.buildDatabase(fineMechanisms, finePoI, fineParam, buildType, paramsGen, offsetGen)

    print '# of further tests: %d' % len(finePoI)
    print 'Calculating distance metrics'
    # metrics = findClosest(testTrace, finePoI, lambda t1,t2: getDistanceVectors(t1,t2)[0])
    metrics = findClosest(testTrace, finePoI, lambda t1,t2: distanceMetric(t1,t2))
    print 'Parameters:'
    print fineParam[metrics[0][1]]
    return finePoI, fineParam, metrics

