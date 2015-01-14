"""Given a set of lengths: crank1, crank2, coupler, base
And also the initial input angle of crank1
We assume the base to be immobile
First calculate the shape of the structure using minimization
Now that we know the joint parameters, we can use forward kinematics to find position of coupler
If necessary, we then offset coupler position to get end effector

Options: 
-h, --help:	Print this text
-d:	Force recalculation of database instead of loading from file
-a: Generate animations of the mechanisms. 
-p: Generates plots of the end effector
-m: Generates plots of the midpoint of the coupler 
"""
import getopt

from math import *
from operator import *
from scipy.optimize import minimize as optimize
import numpy as np
from pylab import *
import matplotlib.pyplot as plt


from Beam2 import Beam
from Line import Line
from data import *
from featureVector import *
from userInput import *
from plotting import *
import construct
import optimizing
import database

mechanismFile = 'mechanisms'
PoIFile = 'PoI'
PoIFileLine = 'PoILine'
parameterLookupFile = 'parameterLookup'
parameterLookupFileLine = 'parameterLookupLine'
optimizationsFolder = 'pics_optimizations/'
NUMPOINTS = database.NUMPOINTS
TRACEMINLEN = 125 #80

mechanisms, PoI, parameterLookup = [], [], []

def start():
	global PoILine
	global parameterLookupLine
	global mechanisms, PoI, parameterLookup

	try: 
		mechanisms = np.load(mechanismFile + '.npy').tolist()
		PoI = np.load(PoIFile + '.npy').tolist()
		parameterLookup = np.load(parameterLookupFile + '.npy').tolist()
	except IOError:
		print '%s %s %s - not found' % (mechanismFile, PoIFile, parameterLookupFile)
		database.buildDatabase(mechanisms, PoI, parameterLookup, 1)
		np.save(mechanismFile, np.array(mechanisms))
		np.save(PoIFile, np.array(PoI))
		np.save(parameterLookupFile, np.array(parameterLookup))
		
	try: 
		PoILine = np.load(PoIFileLine + '.npy')
		parameterLookupLine = np.load(parameterLookupFileLine + '.npy')
	except IOError:
		print '%s %s - not found' % (PoIFileLine, parameterLookupFileLine)
		PoILine, parameterLookupLine = testPoILine()

	# parse command line options
 	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrdpam", ["help"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)
	# process options
	for o, a in opts:
		if o in ("-h", "--help"):
			print __doc__
			sys.exit(0)
		if o in ("-d"):
			mechanisms, PoI, parameterLookup = [], [], []
			database.buildDatabase(mechanisms, PoI, parameterLookup, 1)#,(1,1,2,2,2,2,2,2,1),(1,1,1,1,1))
			np.save(mechanismFile, np.array(mechanisms))
			np.save(PoIFile, np.array(PoI))
			np.save(parameterLookupFile, np.array(parameterLookup))
		if o in ("-p"):
			plotPoI(PoI)
		if o in ("-a"):
			animate(mechanisms)
		if o in ("-m"):
			plotResults_mid(mechanisms)
	# process arguments
	for arg in args:
		process(arg) # process() is defined elsewhere



def calcEndpoint(start, angle, length):
	return (start[0] + length * cos(angle), start[1] + length * sin(angle))


def printFeatureVectors():
	print "Feature Vectors:"
	for mids in [getMids(trace) for trace in PoI]:
		print getFeatureVector(mids)


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
		position = construct.buildState(inCrank, rocker, outCrank, base, a)
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





line, = (None,)
testTrace = inputTest('input.txt')
PoILine = None
parameterLookupLine = None # index of an image is the same as its filename



def demo(test, referenceDatabase, parameterDatabase, metric=getDistanceMetric):
	print 'Finding closest: coarse'
	coarse = optimizing.findClosest(test, referenceDatabase, metric)
	closestCoarse = coarse[0][1]  # 709
	print 'Optimizing on index: %d' % closestCoarse
	optimized = optimizing.optimizeParameters(testTrace, closestCoarse, parameterDatabase, 1, getDistanceMetric)
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
			state = construct.buildStateParam(param, angle=i*angleIncrement)
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
	state = construct.buildStateParam(t_p, angle=1.5)
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
			state = construct.buildStateParamLine(param, angle=i*angleIncrement)
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
	state = construct.buildStateParamLine(param,angle)
	end = 3
	for i in range(3):
		beam = state[i]
		plt.plot((beam[0][0], beam[1][0]), (beam[0][1], beam[1][1]))
	plt.show()
# plt.show()
# return coarse, optimized



start()
