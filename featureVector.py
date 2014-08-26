import numpy as np

def getEdges(mids):
	eges = []
	for index in range(len(mids)):
		edges += [np.subtract(mids[index], mids[index-1])]
	return edges

def getLength(edges):
	return np.sum([sqrt(edge[0]**2 + edge[1]**2 + edge[1]**2) for edge in edges])

def getArea(edges):
	area = 0.0
	for index in range(len(mids)):
		area += np.cross(edges[index], edges[index-1])
	return area

def Distance(mids):
	avgX = np.average(mids[:,0])
	avgY = np.average(mids[:,1])
	avgZ = np.average(mids[:,2])
	cm = [avgX, avgY, avgZ]
	return np.subtract([0,0,0], cm)

def getOrientation(distance, vmax):
	return np.arcsin(np.cross(distance/np.linalg.norm(distance), vmax))