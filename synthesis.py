"""Given a set of lengths: crank1, crank2, coupler, base
And also the initial input angle of crank1
We assume the base to be immobile
First calculate the shape of the structure using minimization
Now that we know the joint parameters, we can use forward kinematics to find position of coupler
If necessary, we then offset coupler position to get end effector

Options: 
-h, --help:	Print this text
-t:	Force recalculation of test data instead of loading from file
-a: Generate animations of the mechanisms. 
-p: Generates plots of the end effector
-m: Generates plots of the midpoint of the coupler 
-c: Generates plots of the principal components
-C: Recomputes the principal components instead of loading from file.
"""

from Beam import Beam
from Line import Line
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
PoIFileLine = 'PoILine'
parameterLookupFile = 'parameterLookup'
parameterLookupFileLine = 'parameterLookupLine'
optimizationsFolder = 'pics_optimizations/'
NUMPOINTS = 128
def start():
	global results
	global components
	global PoI
	global parameterLookup
	global PoILine
	global parameterLookupLine
	try:
		results = loadResults()
	except IOError:
		print 'loadResults no file found'
		results = test()
		printResults(results)
		components = getComponents(results)
		PoI, parameterLookup = testPoI()

	try: 
		components = np.load(componentsFile + '.npy')
	except IOError:
		components = getComponents(results)

	try: 
		PoI = np.load(PoIFile + '.npy')
		parameterLookup = np.load(parameterLookupFile + '.npy')
	except IOError:
		print '%s %s - not found' % (PoIFile, parameterLookupFile)
		PoI, parameterLookup = testPoI()

	try: 
		PoILine = np.load(PoIFileLine + '.npy')
		parameterLookupLine = np.load(parameterLookupFileLine + '.npy')
	except IOError:
		print '%s %s - not found' % (PoIFileLine, parameterLookupFileLine)
		PoILine, parameterLookupLine = testPoILine()

	# parse command line options
 	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrtpamCc", ["help"])
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
			PoI, parameterLookup = testPoI()
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


def test():
	"""Iterate through beam lengths and angles
	beam lengths range from 1 to 4
	angle increments of pi/8 

	results is a list of the start-end coordinates of every beam for each state.
	"""
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
					for angle in range(1,int(NUMPOINTS)):
						a = angle*pi/2.0/(NUMPOINTS/4)
						position = buildState(inCrank, rocker, outCrank, base, a)
						if position:
							r += [position]

					results += [r]
			progress += 1.0/16
			print "%d%%" % (progress*100)
	return results

def testPoI():
	"""Same as test(), but only calculates the PoI not the state
	"""
	results = []
	parameterLookup = []
	progress = 0.0
	numPoI = 25 # Will round down to nearest integer that is double a perfect square
	PoIFineness = trunc(sqrt(numPoI/2.0))
	parameters = []

	print 'Calculating PoI'
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
					for angle in range(1,int(NUMPOINTS)):
						a = angle*pi/2.0/(NUMPOINTS/4)
						position = buildState(inCrank, rocker, outCrank, base, a)
						if position:
							# Iterate through positions of Point of Interest on rocker
							count1 = -1
							for PoIOffset in [x/PoIFineness for x in range(-int(PoIFineness), int(PoIFineness+1))]:
								count1 += 1
								count2 = -1
								for PoIDistance in range(-int(PoIFineness), int(PoIFineness+1)):
									parameters += [[inCrankLength, rockerLength, outCrankLength, baseLength, PoIOffset, PoIDistance]]
									count2 += 1
									rocker.PoIOffset = PoIOffset
									rocker.PoIDistance = PoIDistance
									r[(int(PoIFineness)*2+1)*count1 + count2] += [rocker.PoI()]
					for trace, p in zip(r, parameters):
						if len(trace) > 80:
							results += [trace]
							parameterLookup += [p]
							
			progress += 1.0/16
			print "%d%%" % (progress*100)

	np.save(PoIFile, np.array(results))
	np.save(parameterLookupFile, np.array(parameterLookup))
	return results, parameterLookup

def testPoILine():
	"""Same as testPoI(), but uses the mechanism with a constraining line
	"""
	results = []
	parameterLookup = []
	progress = 0.0
	numPoI = 25 # Will round down to nearest integer that is double a perfect square
	PoIFineness = trunc(sqrt(numPoI/2.0))
	parameters = []

	print 'Calculating PoILine'
	# Iterate through linkage lengths
	for inCrankLength in range(1,5):
		for rockerLength in range(1,5):
			for outCrankLength in range(1,5):
				lineTrack = Line(0,2,3)
				for A in range(1,4):
					A = A+0.0
					lineTrack.a = A
					inCrank = Beam(inCrankLength)
					rocker = Beam(rockerLength)
					outCrank = Beam(outCrankLength)
					r = [[] for _ in range(int(PoIFineness*2+1)**2)]
					for angle in range(1,int(NUMPOINTS)):
						a = angle*pi/2.0/(NUMPOINTS/4)
						position = buildStateLine(inCrank, rocker, outCrank, a, lineTrack)
						if position:
							# Iterate through positions of Point of Interest on rocker
							count1 = -1
							for PoIOffset in [x/PoIFineness for x in range(-int(PoIFineness), int(PoIFineness+1))]:
								count1 += 1
								count2 = -1
								for PoIDistance in range(-int(PoIFineness), int(PoIFineness+1)):
									la, lb, lc = lineTrack.a, lineTrack.b, lineTrack.c
									parameters += [[inCrankLength, rockerLength, outCrankLength, PoIOffset, PoIDistance, A, lb, lc]]
									count2 += 1
									rocker.PoIOffset = PoIOffset
									rocker.PoIDistance = PoIDistance
									r[(int(PoIFineness)*2+1)*count1 + count2] += [rocker.PoI()]
					for trace, p in zip(r, parameters):
						if len(trace) > 100:
							results += [trace]
							parameterLookup += [p]
						
			progress += 1.0/16
			print "%d%%" % (progress*100)

	np.save(PoIFileLine, np.array(results))
	np.save(parameterLookupFileLine, np.array(parameterLookup))
	return results, parameterLookup

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
	Otherwise, returns a list of the start and end xyz-coordinates of  
	each beam in the state.
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

	solveState(0)
	o = optimize(constraint,(0.0,), bounds=((0,2*pi),))
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


def buildStateLine(beam1, coupler, beam2, angle, lineTrack):
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
	coupler.position = [beam1.A*cos(beam1.rotation[0]), beam1.A*sin(beam1.rotation[0]), 0.0]

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

	solveState(0)
	o = optimize(constraint,(0.0,), bounds=((0,2*pi),))
	if o.success and constraint(o.x[0]) < constraintBound:
		solveState(o.x[0])
		# ensure that consistent solution to line intersection is chosen
		if beam2.start()[0] > beam2.end()[0]:
			dist = beam2.A
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

def findClosest(traces, distanceMetric = getDistanceMetric):
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
		inCrank = Beam(p[0])
		rocker = Beam(p[1])
		outCrank = Beam(p[2])
		base = Beam(p[3])
		PoIOffset = p[4]
		PoIDistance = p[5]
		rocker.PoIOffset = PoIOffset
		rocker.PoIDistance = PoIDistance
		for angle in range(1,int(NUMPOINTS)):
			a = angle*pi/2.0/(NUMPOINTS/4)
			position = buildState(inCrank, rocker, outCrank, base, a)
			if not position:
				continue
			trace += [rocker.PoI()]
		return getDistanceMetric(testTrace, trace)

	return optimize(optimizingFunction, p, method = 'Nelder-Mead')


def optimizeParameters(testTrace, index):
	"""	Takes the testTrace and index of a certain trace's parameters in parameterLookup and further optimizes those parameters
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
	delta = 2.0
	numPoI = 25 # Will round down to nearest integer that is double a perfect square
	PoIFineness = trunc(sqrt(numPoI/2.0))
	optimizedParameters = []
	results = []
	progress = 0.0
	parameters = []

	print 'Further Optimizing'
	# Iterate through linkage lengths
	fineness = 5
	print 'Generating finer traces'
	for inCrankLength in range(1,fineness):
		inCrankLength = delta/fineness*(inCrankLength - fineness/2.0)+p[0]
		for rockerLength in range(1,fineness):
			rockerLength = delta/fineness*(rockerLength - fineness/2.0)+p[1]
			for outCrankLength in range(1,fineness):
				outCrankLength = delta/fineness*(outCrankLength - fineness/2.0)+p[2]
				for baseLength in range(1,fineness):
					baseLength = delta/fineness*(baseLength - fineness/2.0)+p[3]

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
					for angle in range(1,int(NUMPOINTS)):
						a = angle*pi/2.0/(NUMPOINTS/4)
						position = buildState(inCrank, rocker, outCrank, base, a)
						if position:
							# Iterate through positions of Point of Interest on rocker
							count1 = -1
							for PoIOffset in [p[4]+x/PoIFineness*delta for x in range(-int(PoIFineness/3), int((PoIFineness+1)/3))]:
								count1 += 1
								count2 = -1
								for PoIDistance in [p[5]+x*delta for x in range(-int(PoIFineness/3), int((PoIFineness+1)/3))]:
									parameters += [[inCrankLength, rockerLength, outCrankLength, baseLength, PoIOffset, PoIDistance]]
									count2 += 1
									rocker.PoIOffset = PoIOffset
									rocker.PoIDistance = PoIDistance
									r[(int(PoIFineness)*2+1)*count1 + count2] += [rocker.PoI()]
					for trace, p in zip(r, parameters):
						if len(trace) > 80:
							results += [trace]
							optimizedParameters += [p]						
			progress += 1.0/16
			print "%d%%" % (progress*100)
	
	print '# of further tests: %d' % len(results)
	print 'Calculating distance metrics'
	metrics = findClosest(results, lambda t1,t2: getDistanceVectors(t1,t2)[0])
	# metrics = findClosest(results, lambda t1,t2: getDistanceMetric(t1,t2))
	print 'Parameters:'
	print optimizedParameters[metrics[0][1]]
	return results, optimizedParameters, metrics


def plotParameters(p):
	""" Used mainly for debugging, as a convenience method to plot using only parameters, so not in plotting.py"""
	trace = traceFromParameters(p)
	plotTrace(trace)

def plotNormalized(trace):
	vmax, vmin = getPrincipalComponents(trace)
	lmax, lmin = getAxisLengths(trace, vmax, vmin)
	trace = normalize(trace, lmax)
	plotTrace(trace)

def traceFromParameters(p):
	trace = []
	inCrank = Beam(p[0])
	rocker = Beam(p[1])
	outCrank = Beam(p[2])
	base = Beam(p[3])
	PoIOffset = p[4]
	PoIDistance = p[5]
	rocker.PoIOffset = PoIOffset
	rocker.PoIDistance = PoIDistance
	for angle in range(1,int(NUMPOINTS)):
		a = angle*pi/2.0/(NUMPOINTS/4)
		position = buildState(inCrank, rocker, outCrank, base, a)
		if not position:
			continue
		trace += [rocker.PoI()]
	return trace

def plotParametersMultiple(ps):
	traces = []
	for p in ps:
		trace = traceFromParameters(p)
		traces += [trace]
	plotTraces(traces)





results = None
components = None
line, = (None,)
PoI = None
parameterLookup = None # index of an image is the same as its filename
testTrace = inputTest('input.txt')
PoILine = None
parameterLookupLine = None # index of an image is the same as its filename

start()


def demo():
	print 'Finding closest: coarse'
	coarse = findClosest(PoI)
	closestCoarse = coarse[0][1]  # 709
	print 'Optimizing on index: %d' % closestCoarse
	optimized = optimizeParameters(testTrace, closestCoarse)
	# showDemo(closestCoarse, optimized)
	return coarse, optimized


def showDemo(closestCoarse, optimized):
	def init():
		for line in lines:
			line.set_data([],[])
		return lines

	def animateTrace(param, frames):
		angleIncrement = 2*pi/frames
		def animate(i):
			i = i+1
			state = buildStateParam(param, angle=i*angleIncrement)
			if not type(state) == list:
				init()
			else:
				mechanismX, mechanismY = [],[]
				for i in range(4):
					beams = state[i]
					mechanismX += [beams[0][0] + offsetDistance[0]]
					mechanismY += [beams[0][1] + offsetDistance[1]]
				lines[0].set_data(mechanismX, mechanismY)
				mechanismX = [x + offsetDistance[0] for x in [state[4][0][0], state[4][1][0],state[5][1][0]]]
				mechanismY = [y + offsetDistance[1] for y in [state[4][0][1], state[4][1][1],state[5][1][1]]]
				lines[1].set_data(mechanismX, mechanismY)
			return tuple(lines)
		return animate

	param = optimized[1][optimized[2][0][1]]
	frames = 25

	f = plt.figure()
	ax1 = plt.subplot(321)
	ax2 = plt.subplot(322, sharex=ax1, sharey=ax1)
	ax3 = plt.subplot(323, sharex=ax1, sharey=ax1)
	ax4 = plt.subplot(324, sharex=ax1, sharey=ax1)

	ax1.scatter([x[0]for x in testTrace], [x[1] for x in testTrace])
	ax1.set_title('Test Trace')

	ax2.scatter([x[0]for x in PoI[closestCoarse-1]], [x[1] for x in PoI[closestCoarse-1]])
	ax2.set_title('Coarse: Closest Trace')
	t = traceFromParameters(param)

	ax3.scatter([x[0]for x in t], [x[1] for x in t])
	ax3.set_title('Optimized')
	t_p = scale(testTrace, param)
	t2 = traceFromParameters(t_p)

	ax4.scatter([x[0]for x in t2], [x[1] for x in t2])
	ax4.set_title('Scaled Optimized')
	state = buildStateParam(t_p, angle=1.5)
	for line in state:
		temp = zip(line[0],line[1])
		ax4.plot(temp[0], temp[1])

	ax5 = plt.subplot(325)
	ax5.scatter([x[0] for x in testTrace], [x[1] for x in testTrace])
	offsetDistance = getDistance(np.array(t2))
	ax5.scatter([x[0]+ offsetDistance[0]for x in t2], [x[1]+ offsetDistance[1] for x in t2])
	ax5.set_title('Animation')
	l, = plt.plot([], [], lw=2,color='black')
	l2, = plt.plot([], [], lw=2,color='red')
	lines = [l, l2]
	plt.xlim(ax1.get_xlim()[0]*2, ax1.get_xlim()[1]*2)
	plt.ylim(ax1.get_ylim()[0]*2, ax1.get_ylim()[1]*2)
	line_ani = animation.FuncAnimation(f, animateTrace(t_p, frames), interval=30, init_func=init, frames=frames)
	mywriter = animation.FFMpegWriter()
	line_ani.save(optimizationsFolder + 'fourbar.mp4', fps=25,writer=mywriter)
	plt.show()
	

def showDemoLine(param):
	def init():
		for line in lines:
			line.set_data([],[])
		return lines

	def animateTrace(param, frames):
		angleIncrement = 2*pi/frames
		def animate(i):
			i = i+1
			state = buildStateParamLine(param, angle=i*angleIncrement)
			if not type(state) == list:
				init()
			else:
				mechanismX, mechanismY = [],[]
				end = 4
				for beam in range(end):
					beams = state[beam]
					mechanismX += [beams[0][0] + offsetDistance[0]]
					mechanismY += [beams[0][1] + offsetDistance[1]]
				lines[0].set_data(mechanismX, mechanismY)
				# mechanismX = [x + offsetDistance[0] for x in [state[end][0][0], state[end][1][0],state[end+1][1][0]]]
				# mechanismY = [y + offsetDistance[1] for y in [state[end][0][1], state[end][1][1],state[end+1][1][1]]]
				# lines[1].set_data(mechanismX, mechanismY)
			return tuple(lines)
		return animate

	frames = 50

	offsetDistance = [0,0,0]

	plt.close()
	fig = plt.figure()
	ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
	l, = plt.plot([], [], 'o-', lw=2,color='black')
	l2, = plt.plot([], [], lw=2,color='red')
	lines = [l, l2]
	line_ani = animation.FuncAnimation(fig, animateTrace(param, frames), interval=30, init_func=init, frames=frames)
	plt.show()

def plotStateLineParam(param, angle =0.1):
	state = buildStateParamLine(param,angle)
	end = 3
	for i in range(3):
		beam = state[i]
		plt.plot((beam[0][0], beam[1][0]), (beam[0][1], beam[1][1]))
	plt.show()
# plt.show()
# return coarse, optimized