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

from Beam import Beam
from math import *
from operator import *
from scipy.optimize import minimize as optimize
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import getopt
from plotting import *
from data import *
from featureVector import *
from userInput import *

componentsFile = 'results_Components'
PoIFile = 'PoI'
def start():
	global results
	global components
	global PoI
	try:
		results = loadResults()
	except IOError:
		results = test()
		printResults(results)
		components = getComponents(results)
		PoI = testPoI()

	try: 
		components = np.load(componentsFile + '.npy')
	except IOError:
		components = getComponents(results)

	try: 
		PoI = np.load(PoIFile + '.npy')
	except IOError:
		PoI = testPoI()

	# parse command line options
 	try:
		opts, args = getopt.getopt(sys.argv[1:], "htpamCcr", ["help"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)
	# process options
	for o, a in opts:
		if o in ("-h", "--help"):
			print __doc__
			sys.exit(0)
		if o in ("-r"):
			results = test()
			printResults(results)
			filterResults(results)
		if o in ("-t"):
			PoI = testPoI()
		if o in ("-p"):
			plotPoI(PoI)
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
	# Iterate through linkage lengths
	for inCrankLength in range(1,5):
		for rockerLength in range(1,5):
			for outCrankLength in range(1,5):
				for baseLength in range(1,5):
					if inCrankLength > baseLength or outCrankLength > baseLength:
						continue
					# Check Grashof condition
					if not (baseLength+rockerLength-inCrankLength-outCrankLength>0 and baseLength-rockerLength-inCrankLength+outCrankLength>0 and -baseLength+rockerLength-inCrankLength+outCrankLength>0):
						continue

					inCrank = Beam(inCrankLength)
					rocker = Beam(rockerLength)
					outCrank = Beam(outCrankLength)
					base = Beam(baseLength)
					r = []
					for angle in range(1,int(numPoints)):
						a = angle*pi/2.0/(numPoints/4)
						position = buildState(inCrank, rocker, outCrank, base, a)
						if position:
							r += [position]

					results += [r]
			progress += 1.0/16
			print "%s%%" % progress
	return results

def testPoI():
	"""Same as above, but only calculates the PoI not the state
	"""
	numPoints = 128.0
	results = []
	progress = 0.0
	numPoI = 25 # Will round down to nearest integer that is double a perfect square
	PoIFineness = trunc(sqrt(numPoI/2.0))

	# Iterate through linkage lengths
	for inCrankLength in range(1,5):
		for rockerLength in range(1,5):
			for outCrankLength in range(1,5):
				for baseLength in range(1,5):
					if inCrankLength > baseLength or outCrankLength > baseLength:
						continue
					# Check Grashof condition
					if not (baseLength+rockerLength-inCrankLength-outCrankLength>0 and baseLength-rockerLength-inCrankLength+outCrankLength>0 and -baseLength+rockerLength-inCrankLength+outCrankLength>0):
						continue

					inCrank = Beam(inCrankLength)
					rocker = Beam(rockerLength)
					outCrank = Beam(outCrankLength)
					base = Beam(baseLength)
					r = [[] for _ in range(int(PoIFineness*2+1)**2)]
					for angle in range(1,int(numPoints)):
						a = angle*pi/2.0/(numPoints/4)
						position = buildState(inCrank, rocker, outCrank, base, a)
						if position:
							# Iterate through positions of Point of Interest on rocker
							count1 = -1
							for PoIOffset in [x/PoIFineness for x in range(-int(PoIFineness), int(PoIFineness+1))]:
								count1 += 1
								count2 = -1
								for PoIDistance in range(-int(PoIFineness), int(PoIFineness+1)):
									count2 += 1
									rocker.PoIOffset = PoIOffset
									rocker.PoIDistance = PoIDistance
									r[(int(PoIFineness)*2+1)*count1 + count2] += [rocker.PoI()]
					for trace in r:
						if len(trace) > 80:
							results += [trace]
							
			progress += 1.0/16
			print "%s%%" % progress
	np.save(PoIFile, np.array(results))
	return results


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

def getComponents(results):
	print "Getting principal components"
	components = []
	for trace in results:
		components += [getPrincipalComponents(getMids(trace))]
	np.save(componentsFile, components)
	return np.array(components)

def printFeatureVectors():
	print "Feature Vectors:"
	for mids in [getMids(trace) for trace in PoI]:
		print getFeatureVector(mids)

results = None
components = None
line, = (None,)
PoI = None

start()

with open('testFeature.txt', 'w') as f:
	div = 32.0
	progress = 0.0
	increment = 0.0
	distances = []
	counter = 0
	for trace in PoI:
		if counter >= increment:
			print "%s%%" % progress
			increment += len(PoI)/div
			progress += 1/div
		counter += 1
		dist = getDistanceMetric(trace, testTrace)
		f.write(str(dist))
		f.write('\n')
		distances += [dist]
	distances = [[distances[index], index+1] for index in range(len(distances))]
	# minIndex = 0
	# for index in range(len(distances)):
	# 	dist = distances[index]
	# 	if dist < distances[minIndex]:
	# 		minIndex = index