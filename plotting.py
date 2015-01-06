# Contains functions to plot results. 

import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from pylab import *
from Beam2 import Beam

# DEPRECATED METHOD
def plotResults_old(results):
	"""Plotting - Assumes results holds the coordinate of the end effector

	results - a list of lists of xyz-coordinates. [[[x1,y1,z1],[x2,y2,z2],...]]
	"""
	print "Plotting results - old method"
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

# DEPRECATED METHOD
def plotResults(results):
	"""Plotting - Assumes results holds the start-end coordinates of each beam.

	results - a list of lists of pairs of xyz-coordinates.
	"""
	folderName = "pics/"
	print "Plotting results: " + folderName
	plt.close()
	color = colorGenerator()
	counter = 0
	for trace in results:
		if len(trace) == 1:
			continue
		ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
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
			plt.savefig(folderName + '%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		plt.close()
		print counter

def plotResults_mid(results):
	"""Plots the midpoint of the coupler
	Assumes results holds the start-end coordinates of each beam.

	results - a list of lists of pairs of xyz-coordinates.
	"""
	folderName = "pics_mid/"
	print "Plotting midpoint of coupler: " + folderName
	counter = 0
	for trace in results:
		ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
		counter += 1
		if len(trace) == 1:
			continue
		for state in trace:
			coupler = state[1]
			plt.scatter((coupler[0][0]+coupler[1][0])/2, (coupler[0][1]+coupler[1][1])/2)
		try:
				plt.savefig(folderName + '/%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		plt.close()
		print counter

def colorGenerator():
	"""Helper function
	Returns a method that cycles through four different colors, returning one each time. For plotting four different beams

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


def plotComponents(results, components):
	folderName = "pics_components/"
	print "Plotting principal components: " + folderName
	plt.close()
	counter = 0
	
	for trace, components in zip(results, components):
		counter += 1
		if len(trace) == 1:
			continue
		points = getMids(trace)
		points = np.array(points)
		
		v1, v2 = components[0], components[1]

		ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
		plt.plot(points[:,0], points[:,1])
		ax.arrow(0,0,v1[0],v1[1], color='red')
		ax.arrow(0,0,v2[0],v2[1])

		try:
			plt.savefig(folderName + '%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		plt.close()
		print counter

# Plots the position of the end-effector
def plotPoI(PoI):
	"""Assumes PoI holds the coordinates of each point.

	PoI - a list of lists of xyz-coordinates.
	"""
	folderName = "pics_PoI/"
	print "Plotting traces of different PoI: " + folderName
	counter = 0
	length = len(PoI)
	checkpoint = length/16.0
	progress = 0
	for trace in PoI:
		plt.clf()
		ax = plt.axes(xlim=(-8,8), ylim=(-8,8))
		counter += 1
		if len(trace) == 1:
			continue
		for point in trace:
			plt.scatter(point[0], point[1])
		try:
				plt.savefig(folderName + '/%d.png' % counter, bbox_inches='tight')
		except AssertionError:
			pass
		if counter >= checkpoint:
			progress += 1/16.0
			checkpoint += length/16.0
			print "%s%%" % progress
	plt.close()


def init():
	line.set_data([], [])
	return line,



# Animates the trace.
def animate(results):
	def animateTrace(results, trace):
		trace = results[trace]
		def animate(i):
			structure = trace[i]
			thisx = [structure[_][0][0] for _ in range(4)]
			thisy = [structure[_][0][1] for _ in range(4)]
			line.set_data(thisx, thisy)
			return line,
		return animate
		
	folderName = "animations/"
	print "Animating: " + folderName
	global line
	for index in range(len(results)):
		frames = len(results[index])
		if frames <= 1:
			continue
		plt.close()
		fig = plt.figure()
		ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
		line, = ax.plot([], [], 'o-', lw=2)
		anim = animation.FuncAnimation(fig, animateTrace(results, index), interval=30, init_func=init, frames=frames)
		anim.save(folderName + str(index) + '.mp4', fps=15)
		print index


def getMids(trace):
	mids = [((state[1][0][0]+state[1][1][0])/2, (state[1][0][1]+state[1][1][1])/2, (state[1][0][2]+state[1][1][2])/2) for state in trace]
	return np.array(mids)


def plotTrace(trace):
	xPoints = [point[0] for point in trace]
	yPoints = [point[1] for point in trace]
	plt.scatter(xPoints, yPoints)
	plt.show()

def plotTraces(traces):
	for trace in traces: 
		xPoints = [point[0] for point in trace]
		yPoints = [point[1] for point in trace]
		plt.scatter(xPoints, yPoints)
	plt.show()