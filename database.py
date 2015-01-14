import math

from Beam2 import Beam
from Line import Line
import construct


NUMPOINTS = 128
TRACEMINLEN = 125

def lengthGenerator(start1=1, end1=4, start2=1, end2=4, start3=1, end3=4, start4=1, end4=4, step=1):
    """Geneartor cycles through tuples that represent a set of length parameters. 
    arguments are the start and end values for each beam, inclusive, and the step to use when iterating.
    """
    assert step>0, 'Step is not positive'

    length1 = start1
    while length1 <= end1:
        length2 = start2
        while length2 <= end2:
            length3 = start3
            while length3 <= end3:
                length4 = start4
                while length4 <= end4:
                    yield (length1, length2, length3, length4)
                    length4 += step
                length3 += step
            length2 += step
        length1 += step

def offsetGenerator(offsetStart=-1, offsetEnd=1, distanceStart=-3, distanceEnd=3, step=1.0):
    assert step>0, 'Step is not positive'

    offset = offsetStart
    while offset <= offsetEnd:
        distance = distanceStart
        while distance <= distanceEnd:
            yield (offset, distance)
            distance += step
        offset += step/3.0



def getBuildFunction(buildType):
    if buildType == 1:
        return construct.buildState
    if buildType == 2:
        return construct.buildStateLine


def isValidParameters(buildType, *args):
    """buildType list:
    1 - 4-bar mechanism
        Expects (length1, length2, length3, length4)

    IS NOT CURRENTLY IMPLEMENTED
    2 - 3-bar mechanism constrained on line
        Expects (length1, length2, length3, lineTrack object)
    """
    if buildType == 1:
        return satisfiesGrashof(*args)
    if buildType == 2:
        return satisfiesLineTrackRadius(*args)


def satisfiesLineTrackRadius(inCrankLength, rockerLength, outCrankLength, lineTrack):
    distance = lineTrack.distanceToPoint(0,0)
    return rockerLength+outCrankLength > inCrankLength+distance


def satisfiesGrashof(inCrankLength, rockerLength, outCrankLength, baseLength):
    """Grashoff Conditions
    baseLength+rockerLength-inCrankLength-outCrankLength > 0 
    baseLength-rockerLength-inCrankLength+outCrankLength > 0
    baseLength+rockerLength-inCrankLength+outCrankLength > 0
    """
    return baseLength+rockerLength-inCrankLength-outCrankLength>0 and baseLength-rockerLength-inCrankLength+outCrankLength>0 and -baseLength+rockerLength-inCrankLength+outCrankLength>0


def countElements(params):
    """params has for (start, end, start, end, ..., step)
    """
    assert len(params)%2 == 1, 'Invalid Parameters'
    index = 0
    count = 1
    while index < len(params)-1:
        temp = params[index]
        tempCount = 0   
        while temp <= params[index+1]:
            tempCount += 1
            temp += params[-1]
        count *= tempCount
        index += 2
    return int(count)


def buildDatabase(mechanismDatabase, traceDatabase, paramDatabase, buildType, paramsGen=(1,4,1,4,1,4,1,4,1), offsetGen=(-1,1,-3,3,1.0)):
    """Builds a database of traces created from 2 fixed cranks.
    traceDatabase - a list to add additional traces to
    paramDatabase - a list to add additional corresponding params to
    paramsGen - a tuple of arguments to pass to lengthGenerator
    offsetGen - a tuple of arguments to pass to offsetGenerator
    buildType - always 1, for now
    """
    startIndex = len(traceDatabase)
    assert startIndex == len(paramDatabase) and startIndex == len(mechanismDatabase), 'mismatched Database indices'

    lengthsList = lengthGenerator(*paramsGen)
    angles = [i*math.pi*2.0/NUMPOINTS for i in range(1, int(NUMPOINTS))]

    buildFunction = getBuildFunction(buildType)
    baseMechanisms = []
    sizeOfLengthList = countElements(paramsGen)
    displayPercent = .05
    displayCounter = displayPercent*sizeOfLengthList
    counter = -1
    for lengths in lengthsList:
        counter += 1
        if counter > displayCounter:
            print '%d%%' % (float(counter)/sizeOfLengthList*100)
            displayCounter += displayPercent*sizeOfLengthList
        inCrankLength, rockerLength, outCrankLength, baseLength = lengths
        if not isValidParameters(buildType, *lengths): continue
        beams = [Beam(length) for length in lengths]
        rocker = beams[1]

        # Calculate if the trace created is 'full' enough
        traceLength = 0
        for angle in angles:
            state = buildFunction(*(beams+[angle]))
            if state:
                traceLength += 1
        if traceLength < TRACEMINLEN:
            continue
        
        traces = []
        mechanisms = []

        # Construct multiple traces simultaneously to avoid redundant buildFunction() calls
        lastAngle = 0.0
        for angle in angles:
            state = buildFunction(*(beams+[angle, lastAngle]))
            angleTraces = []
            angleMechanisms = []
            if state:
                for PoIOffset, PoIDistance in offsetGenerator(*offsetGen):
                    rocker.PoIOffset = PoIOffset
                    rocker.PoIDistance = PoIDistance
                    angleTraces.append(rocker.PoI())
                    angleMechanisms.append(state)
            else:
                for PoIOffset, PoIDistance in offsetGenerator(*offsetGen):
                    angleTraces.append(None)
                    angleMechanisms.append(None)
            traces.append(angleTraces)
            mechanisms.append(angleMechanisms)
            lastAngle = beams[1].rotation[0]

        testLength = len(traces[0])
        for trace in traces:
            assert len(trace) == testLength, 'DEBUG'
        traces = zip(*traces)
        traces = [filter(None, trace) for trace in traces]

        # Add traces and parameters to both databases
        databaseIndex = 0
        for offsetParams in offsetGenerator(*offsetGen):
            traceDatabase.append(traces[databaseIndex])
            mechanismDatabase.append(mechanisms[databaseIndex])
            databaseIndex += 1
            paramDatabase.append(lengths + offsetParams)


def buildDatabaseLine(mechanismDatabase, traceDatabase, paramDatabase, buildType, paramsGen=(1,4,1,4,1,4,1), offsetGen=(-1,1,-3,3,1.0)):
    """Builds a database of traces created from 2 fixed cranks.
    traceDatabase - a list to add additional traces to
    paramDatabase - a list to add additional corresponding params to
    paramsGen - a tuple of arguments to pass to lengthGenerator
    offsetGen - a tuple of arguments to pass to offsetGenerator
    buildType - always 2, for now
    """
    startIndex = len(traceDatabase)
    assert startIndex == len(paramDatabase) and startIndex == len(mechanismDatabase), 'mismatched Database indices'

    lengthsList = lengthGeneratorLine(*paramsGen)
    angles = [i*math.pi*2.0/NUMPOINTS for i in range(1, int(NUMPOINTS))]

    buildFunction = getBuildFunction(buildType)
    baseMechanisms = []
    sizeOfLengthList = 9*countElements(paramsGen)
    displayPercent = .05
    displayCounter = displayPercent*sizeOfLengthList
    counter = -1
    lines = lineGen()

    for slope, intercept in lines:
        if counter > displayCounter:
            print '%d%%' % (float(counter)/sizeOfLengthList*100)
            displayCounter += displayPercent*sizeOfLengthList
        line = (Line(slope, 2, intercept),)
        for lengths in lengthsList:
            counter += 1
            inCrankLength, rockerLength, outCrankLength = lengths
            if not isValidParameters(buildType, *(lengths+line)): continue
            beams = [Beam(length) for length in lengths]
            rocker = beams[1]

            # Calculate if the trace created is 'full' enough
            traceLength = 0
            for angle in angles:
                state = buildFunction(*(beams+[angle,line[0]]))
                if state:
                    traceLength += 1
            if traceLength < TRACEMINLEN:
                continue
            
            traces = []
            mechanisms = []

            # Construct multiple traces simultaneously to avoid redundant buildFunction() calls
            lastAngle = 0.0
            for angle in angles:
                state = buildFunction(*(beams+[angle,line[0], lastAngle]))
                angleTraces = []
                angleMechanisms = []
                if state:
                    for PoIOffset, PoIDistance in offsetGenerator(*offsetGen):
                        rocker.PoIOffset = PoIOffset
                        rocker.PoIDistance = PoIDistance
                        angleTraces.append(rocker.PoI())
                        angleMechanisms.append(state)
                else:
                    for PoIOffset, PoIDistance in offsetGenerator(*offsetGen):
                        angleTraces.append(None)
                        angleMechanisms.append(None)
                traces.append(angleTraces)
                mechanisms.append(angleMechanisms)
                lastAngle = beams[1].rotation[0]

            testLength = len(traces[0])
            for trace in traces:
                assert len(trace) == testLength, 'DEBUG'
            traces = zip(*traces)
            traces = [filter(None, trace) for trace in traces]

            # Add traces and parameters to both databases
            databaseIndex = 0
            for offsetParams in offsetGenerator(*offsetGen):
                traceDatabase.append(traces[databaseIndex])
                mechanismDatabase.append(mechanisms[databaseIndex])
                databaseIndex += 1
                paramDatabase.append(lengths + offsetParams + (slope, intercept))

def lengthGeneratorLine(start1=2, end1=3, start2=1, end2=4, start3=1, end3=4, step=1):
    """Geneartor cycles through tuples that represent a set of length parameters. 
    arguments are the start and end values for each beam, inclusive, and the step to use when iterating.
    """
    assert step>0, 'Step is not positive'

    length1 = start1
    while length1 <= end1:
        length2 = start2
        while length2 <= end2:
            length3 = start3
            while length3 <= end3:
                yield (length1, length2, length3)
                length3 += step
            length2 += step
        length1 += step

def lineGen(start1=0, end1=2, start2=0, end2=2, step=1):
    assert step>0, 'Step is not positive'

    length1 = start1
    while length1 <= end1:
        length2 = start2
        while length2 <= end2:
            yield (length1, length2)
            length2 += step
        length1 += step