from robot import Link
from operator import *
from math import *

class Beam(Link):
	"""Represents a beam, a specific type of link that is a rigid, straight, and has no twist. 
	Calculates positions of endpoints for calculating pin connections
	Position of beam is defined as the startPin
	Travel the length of the Beam to get to the endPin
	The point of interest is the point that we use to trace out a curve. This is measured the perpendicular distnace from the point to the beam, and how far away from the start that intersection is from start. 
	PoIOffset	= distance along the beam, as a ratio of the beam length. Positive offset means traveling towards the endEffector, negative means travelling away. 
	PoIDistance = Perpendicular distance the the PoI. Positive means a clockwise motion, negative means counterclockwise. 
	"""
	zeroThreshold = .00001
	
	
	def __init__(self, length, PoIOffset = 0.0, PoIDistance = 0.0):
		Link.__init__(self, 0,length,0,0)
		self.position = [0.0,0.0,0.0]
		# rotation about Z axis, X axis, Z axis
		self.rotation = [0.0,0.0,0.0]
		# a unit vector describing the direction of the axis that this beam rotates around
		self.axis = [0.0, 0.0, 1.0]
		self.PoIOffset = PoIOffset
		self.PoIDistance = PoIDistance

	def __setattr__(self, name, value):
		if name in Link.fields:
			Link.__setattr__(self, name, value)
		else:
			self.__dict__[name] = value

	def start(self):
		return self.position

	def end(self):
		return travel(self.position, self.rotation[0], self.A)

	def where(self):
		array = [a for a in self.position]
		for i in range(len(array)):
			if abs(array[i]) < self.zeroThreshold:
				array[i] = 0.0
		return array

	def PoI(self):
		intersect = travel(self.position, self.rotation[0], self.PoIOffset*self.A)
		PoI = travel(intersect, self.rotation[0]+pi/2, self.PoIDistance)
		return PoI

	def offsetBeam(self):
		intersect = travel(self.position, self.rotation[0], self.PoIOffset*self.A)
		end = self.PoI()
		return [[list(self.position), list(intersect)], [list(intersect), list(end)]]

def travel(startPos, angle, distance):
	"""Utility function to find relative positions
	Angle is measured from horizontal pointing right."""
	ret = map(add, startPos, [distance*x for x in (cos(angle), sin(angle), 0.0)])
	return ret

