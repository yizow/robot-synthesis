# DEPRECATED CLASS NO LONGER USED
import numpy as np
from matplotlib.mlab import PCA

class Curve():
	def __init__(points):
		mid = np.array(points)
		if mid.shape[0] < mid.shape[1]:
			mid = mid.T
		pca = PCA(mid)

		X = pca.Y[:,0]
		Y = pca.Y[:,1]
		origX = pca.a[:,0]
		origY = pca.a[:,1]
		components = (pca.Wt[0], pca.Wt[1])
		