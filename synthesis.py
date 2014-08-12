"""Given a set of lengths: crank1, crank2, coupler, base
And also the initial input angle of crank1
We assume the base to be immobile
First calculate the shape of the structure using minimization
Now that we know the joint parameters, we can use forward kinematics to find position of coupler
If necessary, we then offset coupler position to get end effector

Options: 
-h, --help:	Print this text
-t:	Force recalculation of test data instead of loading from file
"""

from math import *
from robot import Link
from operator import *
from scipy.optimize import minimize as optimize
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.path import Path
import getopt

def start():
	global results
	_p = False;

	# parse command line options
 	try:
		opts, args = getopt.getopt(sys.argv[1:], "htp", ["help"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)
	# process options
	for o, a in opts:
		if o in ("-h", "--help"):
			print __doc__
			sys.exit(0)
		if o in ("-t"):
			results = test()
			printResults(results)
		if o in ("-p"):
			_p = True
	# process arguments
	for arg in args:
		process(arg) # process() is defined elsewhere


	try:
		results = loadResults()
	except IOError:
		results = test()
		printResults(results)

	if _p:
		plotResults(results)

def calcEndpoint(start, angle, length):
	return (start[0] + length * cos(angle), start[1] + length * sin(angle))


def euler(theta):
	"""theta is rotation about z-axis. Only allowed in 2-D
	"""
	return tr2eul(r2t(rotz(theta)))


class Beam(Link):
	"""Represents a beam, a specific type of link that is a rigid, straight, and has no twist. 
	Calculates positions of endpoints for calculating pin connections
	Position of beam is defined as the startPin
	Travel the length of the Beam to get to the endPin
	"""
	zeroThreshold = .00001
	
	
	def __init__(self, length, endEffector = [0.0,0.0,0.0]):
		"""endEffector is vector pointing from start to endEffector position, if self is at 0 degrees rotation
		"""
		Link.__init__(self, 0,length,0,0)
		self.position = [0.0,0.0,0.0]
		# rotation about Z axis, X axis, Z axis
		self.rotation = [0.0,0.0,0.0]
		# a unit vector describing the direction of the axis that this beam rotates around
		self.axis = [0.0, 0.0, 1.0]
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


def pinConnection(beam1, beam2):
	"""Uses the end pin on beam1 and the start pin of beam2
	"""
	posConstraint = map(abs, map(sub, beam1.end(), beam2.start()))
	# rotConstraint = map(abs, map(sub, beam1.position, beam2.position))
	rotConstraint = cross(beam1.axis, beam2.axis)
	return posConstraint + rotConstraint


def buildConstraint(beam1, coupler, beam2, base):
	return np.add(np.add(pinConnection(beam1, coupler), pinConnection(coupler, beam2)), np.add(pinConnection(beam2, base), pinConnection(base, beam1)))

def buildState(beam1, coupler, beam2, base, angle):
	"""Returns none if optimization fails to find a zero (impossible state)
	Otherwise, returns the xyz-coordinates of the: 
	beam2
	"""
	constraintBound = 0.001

	base.position = [base.A,0.0,0.0]
	base.rotation = [pi,0.0,0.0]
	beam1.position = [0.0,0.0,0.0]
	beam1.rotation = [angle,0.0,0.0]
	coupler.position = [beam1.A*cos(beam1.rotation[0]), beam1.A*sin(beam1.rotation[0]), 0.0]

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
	if o.success and constraint(o.x[0]) < constraintBound:
		solveState(o.x[0])
		# print constraint(o.x[0])
		return [[list(beam.start()), list(beam.end())] for beam in (beam1, coupler, beam2, base)]
	else:
		return None


# b1 = Beam(3)
# c = Beam(3)
# b2 = Beam(3)
# b = Beam(3)
# a = pi/6

def test():
	"""Iterate through beam lengths and angles
	beam lengths range from 1 to 4
	angle increments of pi/8 

	results is a list of the start-end coordinates of every beam for each state.
	"""
	angIncrement = 128.0
	results = []
	progress = 0.0
	for b1Length in range(1,5):
		b1 = Beam(b1Length)
		for couplerLength in range(1,5):
			c = Beam(couplerLength)
			for b2Length in range(1,5):
				b2 = Beam(b2Length)
				for baseLength in range(1,5):
					b = Beam(baseLength)
					r = []
					for angle in range(1,int(angIncrement)):
						a = angle*pi/2.0/(angIncrement/4)
						position = buildState(b1,c,b2,b,a)
						if position:
							r += [position]

					results += [r]
			progress += 1.0/16
			print "%s%%" % progress
	return results


def printResults(results):
	with open('results.txt', 'w') as f:
		for items in results:
			s = ""
			for coordinate in items:
				s += "%s:" % coordinate
			s = s[:-1] + "\n" # drop the last colon and new line
			f.write(s)


def loadResults():
	with open('results.txt', 'r') as f:
		results = []
		
		while True:
			line = f.readline()
			
			# File End. Return
			if line == '':	
				return results

			trace = []
			for coordinate in line.split(':'):
				# Remove newline character, brackets, and commas; then split each coordinate into values
				coordinate = coordinate[:-1].translate(None, '[],').split()
				# Convert to floats and transform back into start-end pairs; then add to trace
				pairs = []
				for i in range(len(coordinate)/6):
					index = i * 6
					try:
						for add in range(6):
							coordinate[index + add] = float(coordinate[index+add])
					except ValueError:
						print coordinate[i]
					startEnd = [[coordinate[index], coordinate[index+1], coordinate[index+2]], [coordinate[index+3],coordinate[index+4],coordinate[index+5]]]
					pairs += [startEnd]

				trace += [pairs]

			results += [trace]




def plotResults_old(results):
	"""Plotting - Assumes results holds the coordinate of the end effector

	results - a list of lists of xyz-coordinates. [[[x1,y1,z1],[x2,y2,z2],...]]
	"""
	for coordinates in results:
		x = []
		y = []
		color = []
		for coordinate in coordinates:
			x += [coordinate[0]]
			y += [coordinate[1]]
			color += ['r']
		fig = plt.figure()
		ax = fig.add_subplot(111)

		scatter(x,y, s=100 ,marker='o', c=color)

		# [ plot( [dot_x,dot_x] ,[0,dot_y]) for dot_x,dot_y in zip(x,y) ] 
		# [ plot( [0,dot_x] ,[dot_y,dot_y]) for dot_x,dot_y in zip(x,y) ]

		left,right = ax.get_xlim()
		low,high = ax.get_ylim()
		arrow( left, 0, right -left, 0)
		arrow( 0, low, 0, high-low) 
		grid()
		show()

def plotResults(results):
	"""Plotting - Assumes results holds the start-end coordinates of each beam.

	results - a list of lists of pairs of xyz-coordinates.
	"""
	color = colorGenerator()
	counter = 0
	for trace in results:
		counter += 1
		for state in trace:
			# xlist = []
			# ylist = []
			# for pair in state:
			# 	xlist.extend((pair[0][0], pair[1][0]))
			# 	ylist.extend((pair[0][1], pair[1][1]))
			# 	# Appending None improves performance
			# 	xlist.append(None)
			# 	ylist.append(None)
			for pair in state:
				plt.plot([pair[0][0], pair[1][0]],[pair[0][1], pair[1][1]], color=color())
				plt.plot(None,None)

		# plt.show()
		try:
			plt.savefig('pics/%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		plt.close()
		print counter

def plotResults_mid(results):
	"""Plots the midpoint of the coupler
	"""
	counter = 0
	for trace in results:
		counter += 1
		if len(trace) == 1:
			continue
		for state in trace:
			coupler = state[1]
			plt.scatter((coupler[0][0]+coupler[1][0])/2, coupler[0][1]+coupler[1][1])
		try:
				plt.savefig('pics_mid/%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		plt.close()
		print counter

def colorGenerator():
	"""Returns a method that cycles through four different colors, returning one each time. For plotting four different beams

	Colors:
	red, blue, green, black
	"""
	colors = ['red', 'blue', 'green', 'black']
	counter = {'index':-1}
	def cycle():
		# nonlocal does not exist in python 2.7
		# nonlocal counter
		counter['index'] = (counter['index'] + 1) %4
		return colors[counter['index']]
	return cycle

results = None
start()
# results = test()
# printResults(results)
# plotResults(results)
# loadResults()