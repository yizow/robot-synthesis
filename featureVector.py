import numpy as np
import math

def getDistanceMetric(trace1, trace2):
	f1 = getFeatureVector(trace1)
	f2 = getFeatureVector(trace2)
	distanceVectors = getDistanceVectors(trace1, trace2)
	diff = np.subtract(f1,f2)
	diff = diff.tolist()
	diff += [distanceVectors[0], distanceVectors[1]]
	distance = np.dot(np.transpose(diff),diff)
	return distance


def getFeatureVector(mids):
	vmax, vmin = getPrincipalComponents(mids)
	lmax, lmin = getAxisLengths(mids, vmax, vmin)
	mids = normalize(mids, lmax)
	edges = getEdges(mids)
	features = []
	features.append(getLength(edges))
	features.append(getArea(mids))
	features.append(lmin/lmax)
	distance = getDistance(mids)
	# Explicit norm calculated for speed
	features.append(math.sqrt(distance[0]**2+distance[1]**2+distance[2]**2))
	features.append(getOrientation(distance, vmax))
	features.append(getNumIntersections(mids))
	return features

def normalize(mids, lmax):
	return np.array([mid/lmax for mid in mids])


def getEdges(mids):
	edges = []
	for index in range(len(mids)):
		edges += [np.subtract(mids[index], mids[index-1])]
	return np.array(edges)

def getLength(edges):
	return np.sum([np.sqrt(edge[0]**2 + edge[1]**2 + edge[1]**2) for edge in edges])

def getArea(mids):
	area = 0.0
	for index in range(len(mids)):
		# Explicit norm calculated for speed
		v = np.cross(mids[index], mids[index-1])
		area += math.sqrt(v[0]**2 + v[1]**2 + v[1]**2)
	return area

def getDistance(mids):
	avgX = np.average(mids[:,0])
	avgY = np.average(mids[:,1])
	avgZ = np.average(mids[:,2])
	cm = [avgX, avgY, avgZ]
	return np.subtract([0,0,0], cm)

def getOrientation(distance, vmax):
	# Explicit norm calculated for speed
	v = np.cross(distance/math.sqrt(distance[0]**2 + distance[1]**2 + distance[1]**2), vmax)
	return np.arcsin(math.sqrt(v[0]**2 + v[1]**2 + v[1]**2))

def getAxisLengths(mids, v1, v2):
	"""Get the length of the curve on the major and minor axis, returned in that order.
	"""
	transform = np.hstack((v1.reshape(3,1), v2.reshape(3,1)))
	transformed = transform.T.dot(np.array(mids).T)
	ranges = np.ptp(transformed, axis=1)
	ranges = [(ranges[0], v1), (ranges[1], v2)]
	if ranges[0][0] > ranges[1][0]:
		ranges = [ranges[1], ranges[0]]
	lmax = ranges[0][0]/2.0
	lmin = ranges[1][0]/2.0
	return (lmax, lmin)

def getPrincipalComponents(mids):
	eVal, eVec = getEig(mids)
	v1, v2 = findPrincipalComponents(eVal, eVec)
	return v1, v2

def getEig(trace):
	trace = np.array(trace)
	meanx = np.average(trace[:,0])
	meany = np.average(trace[:,1])	
	meanz = np.average(trace[:,2])
	correctedX = [value-meanx for value in (trace[:,0])] 
	correctedY = [value-meany for value in (trace[:,1])] 
	correctedZ = [value-meanz for value in (trace[:,2])] 

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

def getNumIntersections(mids):
	length = len(mids)
	counter = 0
	for index in range(length-1):
		A, B = mids[index], mids[index+1]
		for increment in range (2, length-1):
			increment -= length
			C, D = mids[index+increment], mids[index+increment+1]
			if intersects(A,B,C,D):
				counter += 1
				# print str(A) + " - " + str(B)
				# print str(C) + " - " + str(D)
				# print
	return counter

def intersects(A, B, C, D):
	"""Checks for intersection between line segments AB and CD
	"""
	if max(A[0],B[0]) < min(C[0], D[0]) or min(A[0],B[0]) > max(C[0], D[0]):
		return False
	# Two lines are represent by f1(x) = A1*x + b1, f2(x) = A2*x + b2, 
	try:
		A1 = (A[1]-B[1])/(A[0]-B[0])
		A2 = (C[1]-D[1])/(C[0]-D[0])
		b1 = A[1]-A1*A[0]
		b2 = C[1]-A2*C[0]
	except ZeroDivisionError:	# We already checked x-bounds, so this means we have two vertical lines. Check y-bounds
		if max(A[1],B[1]) < min(C[1], D[1]) or min(A[1],B[1]) > max(C[1], D[1]):
			return False
		else: 
			return True

	if A1 == A2:	# Parallel Lines
		return False
	X = (b2-b1)/(A1-A2)

	if X < max(min(A[0],B[0]), min(C[0],D[0])) or X > min(max(A[0],B[0]), max(C[0],D[0])):
	   return False; # intersection is out of bound
	else:
 		return True;

def getDistanceVectors(T1, T2):
	# ensure len(T1) < len(T2)
	if len(T1) > len(T2):
		T1, T2 = T2, T1

	T1Curvature, T2Curvature = discreteCurvature(T1), discreteCurvature(T2)
	dPos, dAngle = [], []
	for increment in range(len(T2)):
		posSum, angleSum = 0.0,0.0
		for index in range(len(T1)):
			difference = np.subtract(T1[index],T2[(index+increment)%len(T2)])
			# explicit norm calculated for speed
			posSum += difference[0]**2+difference[1]**2+difference[2]**2
			angleDifference = T1Curvature[index]-T2Curvature[(index+increment)%len(T2)]
			angleSum += angleDifference[0]**2+angleDifference[1]**2+angleDifference[2]**2
		dPos += [posSum]
		dAngle += [angleSum]

	dPos = math.sqrt(min(dPos)/len(T1))
	dAngle = math.sqrt(min(dAngle))
	return dPos, dAngle

def discreteCurvature(curve):
	edges = getEdges(curve)
	# Explicit norm calculated for speed
	normEdges = [edge/math.sqrt(edge[0]**2+edge[1]**2+edge[2]**2) for edge in edges]
	curvatures = []
	for index in range(len(edges)):
		curvatures += [np.cross(2*normEdges[index-1], normEdges[index])/(1+np.dot(normEdges[index-1], normEdges[index]))]
	return np.array(curvatures)