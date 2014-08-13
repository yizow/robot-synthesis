import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(-5,5), ylim=(-5,5))
line, = ax.plot([], [], 'o-', lw=2)

def init():
	line.set_data([], [])
	return line,

def animate(results, trace):
	trace = results[trace]
	def animate(i):
		structure = trace[i]
		thisx = [structure[_][0][0] for _ in range(4)]
		thisy = [structure[_][0][1] for _ in range(4)]
		line.set_data(thisx, thisy)
		return line,
	return animate