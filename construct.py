from math import *
from scipy.optimize import minimize as optimize
import numpy as np
from operator import *
from Beam2 import Beam

def pinConnection(beam1, beam2):
    """Uses the end pin on beam1 and the start pin of beam2
    """
    posConstraint = map(abs, map(sub, beam1.end(), beam2.start()))
    # rotConstraint = map(abs, map(sub, beam1.position, beam2.position))
    rotConstraint = np.linalg.norm(np.cross(beam1.axis, beam2.axis))
    return posConstraint + [rotConstraint]


def buildConstraint(beam1, coupler, beam2, base):
    return np.add(np.add(pinConnection(beam1, coupler), pinConnection(coupler, beam2)), np.add(pinConnection(beam2, base), pinConnection(base, beam1)))

def buildState(beam1, coupler, beam2, base, angle, lastAngle=0.0):
    """Returns none if optimization fails to find a zero (impossible state)
    Otherwise, returns a list of the start and end xyz-coordinates of  
    each beam in the state.
    """
    constraintBound = 0.001

    base.position = [base.length,0.0,0.0]
    base.rotation = [pi,0.0,0.0]
    beam1.position = [0.0,0.0,0.0]
    beam1.rotation = [angle,0.0,0.0]
    coupler.position = [beam1.length*cos(beam1.rotation[0]), beam1.length*sin(beam1.rotation[0]), 0.0]

    def constraint(theta):
        """beam1 and base state are known. 
        Solve the optimization problem to find states of coupler and beam2
        Input is z-rotation of coupler
        """
        solveState(theta)
        c = buildConstraint(beam1, coupler, beam2, base)
        # print [b.position for b in (beam1, coupler, beam2, base)]
        ret = np.dot(np.transpose(c), c)
        # print ret
        return ret

    def solveState(theta):
        coupler.rotation = [theta,0.0,0.0]
        beam2.position = [a for a in coupler.end()]
        dx = beam2.position[0]-base.position[0]
        dy = beam2.position[1]-base.position[1]
        hyp = sqrt(dx**2 + dy**2)
        beam2.rotation = [acos(dx/hyp)-pi,0.0,0.0]

    solveState(lastAngle)
    o = optimize(constraint,(lastAngle,))
    if o.success and constraint(o.x[0]) < constraintBound:
        solveState(o.x[0])
        # print constraint(o.x[0])
        ret = [[list(beam.start()), list(beam.end())] for beam in (beam1, coupler, beam2, base)]
        ret += coupler.offsetBeam()
        return ret
    else:
        return None

def buildStateParam(param, angle = 0.1):
    """p = [inCrankLength rockerLength outCrankLength baseLength PoIOffset PoIDistance]"""
    beam1 = Beam(param[0])
    coupler = Beam(param[1])
    coupler.PoIOffset = param[4]
    coupler.PoIDistance = param[5]
    beam2 = Beam(param[2])
    base = Beam(param[3])
    return buildState(beam1, coupler, beam2, base, angle)

def buildStateParamLine(param, angle = 0.1):
    """p = [inCrankLength rockerLength outCrankLength PoIOffset PoIDistance a b c]"""
    beam1 = Beam(param[0])
    coupler = Beam(param[1])
    coupler.PoIOffset = param[3]
    coupler.PoIDistance = param[4]
    beam2 = Beam(param[2])
    lineTrack = Line(param[5], param[6], param[7])
    return buildStateLine(beam1, coupler, beam2, angle, lineTrack)


def buildStateLine(beam1, coupler, beam2, angle, lineTrack, lastAngle=0.0):
    """This state has the outcrank sliding on a line, 
    instead of a fixed point.
    lineTrack is a Line object that represents the track the
    outcrank is constrained on.

    Returns none if optimization fails to find a zero (impossible state)
    Otherwise, returns a list of the start and end xyz-coordinates of  
    each beam in the state.
    """
    constraintBound = 0.1

    beam1.position = [0.0,0.0,0.0]
    beam1.rotation = [angle,0.0,0.0]
    coupler.position = [beam1.length*cos(beam1.rotation[0]), beam1.length*sin(beam1.rotation[0]), 0.0]

    def buildConstraint(beam1, coupler, beam2):
        return np.append(np.add(pinConnection(beam1, coupler), pinConnection(coupler,beam2)), lineTrack.distanceToPoint(beam2.end()[0],beam2.end()[1]))

    def constraint(theta):
        """beam1 and base state are known. 
        Solve the optimization problem to find states of coupler and beam2
        Input is z-rotation of coupler
        """
        solveState(theta)
        c = buildConstraint(beam1, coupler, beam2)
        ret = np.dot(np.transpose(c), c)
        return ret

    def solveState(theta):
        coupler.rotation = [theta,0.0,0.0]
        beam2.position = [a for a in coupler.end()]
        dx = beam2.start()[0]-beam2.end()[0]
        dy = beam2.start()[1]-beam2.end()[1]
        hyp = sqrt(dx**2 + dy**2)
        beam2.rotation = [acos(dx/hyp)-pi,0.0,0.0]

    solveState(lastAngle)
    o = optimize(constraint,(lastAngle,))
    if o.success and constraint(o.x[0]) < constraintBound:
        solveState(o.x[0])
        # ensure that consistent solution to line intersection is chosen
        if beam2.start()[0] > beam2.end()[0]:
            dist = beam2.length
            newX = beam2.start()[0] + dist/math.sqrt(1+lineTrack.slope()**2)
            newY = lineTrack.solveForY(newX)
            dx = newX - beam2.start()[0]
            dy = newY - beam2.start()[1]
            hyp = sqrt(dx**2 + dy**2)
            beam2.rotation = [acos(dx/hyp)-pi,0.0,0.0]
        ret = [[list(beam.start()), list(beam.end())] for beam in (beam1, coupler, beam2)]
        ret += coupler.offsetBeam()
        return ret
    else:
        return None