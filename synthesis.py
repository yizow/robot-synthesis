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
import getopt
from plotting import *
from data import *

componentsFile = 'results_Components'
def start():
	global results
	global components
	try:
		results = loadResults()
	except IOError:
		results = test()
		printResults(results)
		components = getComponents(results)

	try: 
		components = np.load(componentsFile + 'npy')
	except IOError:
		components = getComponents(results)

	# parse command line options
 	try:
		opts, args = getopt.getopt(sys.argv[1:], "htpamCc", ["help"])
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
			plotResults(results)
		if o in ("-a"):
			animate(results)
		if o in ("-m"):
			plotResults_mid(results)
		if o in ("-C"):
			components = getComponents(results)
		if o in ("-c"):
			plotComponents(results, components)
	# process arguments
	for arg in args:
		process(arg) # process() is defined elsewhere



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
	rotConstraint = np.linalg.norm(cross(beam1.axis, beam2.axis))
	return posConstraint + [rotConstraint]


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
	o = optimize(constraint,(0.0,), bounds=((0,2*pi),))
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
	numPoints = 128.0
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
					for angle in range(1,int(numPoints)):
						a = angle*pi/2.0/(numPoints/4)
						position = buildState(b1,c,b2,b,a)
						if position:
							r += [position]

					results += [r]
			progress += 1.0/16
			print "%s%%" % progress
	return results


def getLengths(trace):
	"""Get the eigenvectors of the major and minor axis, returned in that order.
		The eigenvectors are scaled by dividing by the length of the trace in the major axis
	"""
	mids = getMids(trace)
	eVal, eVec = getEig(mids)
	v1, v2 = findPrincipalComponents(eVal, eVec)
	transform = np.hstack((v1.reshape(3,1), v2.reshape(3,1)))
	transformed = transform.T.dot(np.array(mids).T)
	ranges = np.ptp(transformed, axis=1)
	ranges = [(ranges[0], v1), (ranges[1], v2)]
	ranges.sort()
	ranges.reverse()
	lmax = ranges[0][0]
	return np.array((ranges[0][1]/lmax, ranges[1][1]/lmax))

def getComponents(results):
	print "Getting principal components"
	components = []
	for trace in results:
		components += [getLengths(trace)]
	np.save(componentsFile, components)
	return np.array(components)


def getEig(points):
	points = np.array(points)
	meanx = np.average(points[:,0])
	meany = np.average(points[:,1])	
	meanz = np.average(points[:,2])
	correctedX = [value-meanx for value in (points[:,0])] 
	correctedY = [value-meany for value in (points[:,1])] 
	correctedZ = [value-meanz for value in (points[:,2])] 

	data = np.array([correctedX, correctedY, correctedZ])
	covData = np.cov(data)
	eigenvalues, eigenvectors = np.linalg.eig(covData)

	return eigenvalues, eigenvectors

def findPrincipalComponents(eigenvalues, eigenvectors):
	"""Given two numpy arrays, one of eigenvalues, the other of the corresponding eigenvalues, 
		This function returns the 2 eigenvectors with the largest eigenvalues
	"""
	value1, value2 = -1.0, -1.0 	#value1 >= value2
	vec1, vec2 = None, None

	# Walk through the eigenvalues and record the two largest 
	for index in range(len(eigenvalues)):
		value = eigenvalues[index]
		if value > value1:
			value2 = value1
			value1 = value
			vec2 = vec1
			vec1 = eigenvectors[index]
		elif value > value2:
			value2 = value
			vec2 = eigenvectors[index]

	return vec1, vec2


results = None
components = None
line, = (None,)

start()

	# plt.show()
# results = test()
# printResults(results)
# plotResults(results)
# loadResults()