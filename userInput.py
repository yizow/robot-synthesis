from numpy import linspace,exp
from numpy.random import randn
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# User defines a number of points. These points will be saved as a sequential list of X coordinate, and Y coordinates
# Compute splines for each of these lists individually.
# Plot the resulting splines as parametric equations to get the user's curve
# Generate the list of xyz points by 


def getPoints(spline, linspace):
	points = spline(linspace)
	points = points[:len(points)]
	return points

def inputTest(filename):
	numTracePoints = 128
	# File expects input as space delimited integers on three lines
	# First line contains X coordinate, second line contains Y, third contains Z
	with open(filename, 'r') as f:
		X = [float(i) for i in f.readline().split()]
		Y = [float(i) for i in f.readline().split()]
		Z = [float(i) for i in f.readline().split()]
		# add back last point to ensure closed path
		X += [X[0]]
		Y += [Y[0]]
		Z += [Z[0]]
	# Check that X,Y,Z have same length. Otherwise, error
	if len(X) != len(Y) or len(Y) != len(Z):
		raise Exception("Input coordinate lists do not have same length")
	numPoints = len(X)
	t = linspace(0,numPoints-1, numTracePoints)
	xs = UnivariateSpline(range(numPoints), X, s=0)
	ys = UnivariateSpline(range(numPoints), Y, s=0)
	zs = UnivariateSpline(range(numPoints), Z, s=0)
	xPoints = getPoints(xs, t)
	yPoints = getPoints(ys, t)
	zPoints = getPoints(zs, t)
	trace = zip(xPoints, yPoints, zPoints)
	return trace
