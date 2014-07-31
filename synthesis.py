from math import *
from robot import Link
from operator import *
from scipy.optimize import minimize as optimize
import numpy as np

def calcEndpoint(start, angle, length):
	return (start[0] + length * cos(angle), start[1] + length * sin(angle))

# theta is rotation about z-axis. Only allowed in 2-D
def euler(theta):
	return tr2eul(r2t(rotz(theta)))

# Given a set of lengths: crank1, crank2, coupler, base
# And also the initial input angle of crank1
# We assume the base to be immobile
# First calculate the shape of the structure using minimization
# Now that we know the joint parameters, we can use forward kinematics to find position of coupler
# If necessary, we then offset coupler position to get end effector

# represents a beam, a specific type of link that is a rigid, straight, and has no twist. 
# Calculates positions of endpoints for calculating pin connections
# Position of beam is defined as the startPin
# Travel the length of the Beam to get to the endPin
class Beam(Link):
	zeroThreshold = .00001
	
	# endEffector is vector pointing from start to endEffector position, if self is at 0 degrees rotation
	def __init__(self, length, endEffector = [0.0,0.0,0.0]):
		Link.__init__(self, 0,length,0,0)
		self.position = [0.0,0.0,0.0]
		# rotation about Z axis, X axis, Z axis
		self.rotation = [0.0,0.0,0.0]
		self.endEffector = [_ for _ in endEffector]

	def __setattr__(self, name, value):
		if name in Link.fields:
			Link.__setattr__(self, name, value)
		else:
			self.__dict__[name] = value

	def start(self):
		return self.position

	def end(self):
		angle = self.rotation[0]
		ret = map(add, self.position, [self.A*x for x in (cos(angle), sin(angle), 0.0)])
		return ret

	def where(self):
		array = [a for a in self.position]
		for i in range(len(array)):
			if abs(array[i]) < self.zeroThreshold:
				array[i] = 0.0
		return array


# Uses the end pin on beam1 and the start pin of beam2
def pinConnection(beam1, beam2):
	posConstraint = map(abs, map(sub, beam1.end(), beam2.start()))
	# rotConstraint = map(abs, map(sub, beam1.position, beam2.position))
	rotConstraint = [0,0,0]
	return posConstraint + rotConstraint

def buildConstraint(beam1, coupler, beam2, base):
	return np.add(np.add(pinConnection(beam1, coupler), pinConnection(coupler, beam2)), np.add(pinConnection(beam2, base), pinConnection(base, beam1)))

# def constraint

def buildState(beam1, coupler, beam2, base, angle):
	base.position = [base.A,0.0,0.0]
	base.rotation = [pi,0.0,0.0]
	beam1.position = [0.0,0.0,0.0]
	beam1.rotation = [angle,0.0,0.0]
	coupler.position = [beam1.A*cos(beam1.rotation[0]), beam1.A*sin(beam1.rotation[0]), 0.0]

	# beam1 and base state are known. 
	# Solve the optimization problem to find states of coupler and beam2
	# Input is z-rotation of coupler
	def constraint(theta):
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

	def findBounds():
		dx = coupler.position[0]-base.position[0]
		dy = coupler.position[1]
		hyp = sqrt(dx**2 + dy**2)
		# if hyp == 0:
		# 	print coupler.position
		# 	print base.position
		# 	print beam1.position
		# 	print beam1.rotation
		angle = acos(dx/hyp)
		# return ((-pi/2,angle),)
		return ((beam1.rotation[0]-pi,angle),)

	solveState(0)
	o = optimize(constraint,(0.0,), bounds=findBounds())
	# print o
	solveState(o.x[0])
	return [_ for _ in beam2.position]


# b1 = Beam(3)
# c = Beam(3)
# b2 = Beam(3)
# b = Beam(3)
# a = pi/6

def test():
	results = []
	# Iterate through beam lengths and angles
	# beam lengths range from 1 to 4
	# angle increments of pi/8 
	for b1Length in range(1,5):
		b1 = Beam(b1Length)
		for couplerLength in range(1,5):
			c = Beam(couplerLength)
			for b2Length in range(1,5):
				b2 = Beam(b2Length)
				for baseLength in range(1,5):
					b = Beam(baseLength)
					r = []
					for angle in range(1,16):
						a = angle*pi/2.0/4.0

						r += [buildState(b1,c,b2,b,a)]

					results += [r]

	with open('results.txt', 'w') as f:
		for item in results:
			f.write("%s\n" % item)

	return results

test()