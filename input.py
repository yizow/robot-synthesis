from numpy import linspace,exp
from numpy.random import randn
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# User defines a number of points. These points will be saved as a sequential list of X coordinate, and Y coordinates
# Compute splines for each of these lists individually.
# Plot the resulting splines as parametric equations to get the user's curve
# Generate the list of xyz points by 

repeat = 4
X = [1,0,-1,0]
X*=repeat
Y = [0,1,0,-1]
Y*=repeat
# Can easily extend this spline representation to 3 dimensions
# Z = [0,0,0,0]
numPoints = len(X)

def getPoints(spline, linspace):
	points = spline(linspace)
	points = points[:len(points)]
	return points

t = linspace(0,numPoints-1, (numPoints-1)*10)
xs = UnivariateSpline(range(numPoints), X, s=0)
ys = UnivariateSpline(range(numPoints), Y, s=0)
# zs = UnivariateSpline(t, Z, s=1)
t = t[:len(t)/repeat]
xPoints = getPoints(xs, t)
yPoints = getPoints(ys, t)

# Plot the user spline for testing
print xPoints
print yPoints
plt.plot(t,xPoints,'r')
plt.plot(t,yPoints,'g')
plt.show()
plt.figure()
plt.scatter(xPoints, yPoints)
plt.show()