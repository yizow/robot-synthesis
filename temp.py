import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_line(num, data, data2, line, line2):
    line.set_data(data[...,:num])
    line2.set_data(data2[...,:num])
    return line,line2 

fig1 = plt.figure()

data = np.array([[1,2,3,4],[2,4,1,2]])/10.0
data2 = np.array([[2,4,1,2],[1,2,3,4]])/10.0
ax1=plt.subplot(311)
ax3=plt.subplot(313)
l, = plt.plot([], [])
l2, = plt.plot([], [], 'r')
ax2=plt.subplot(312)

line_ani = animation.FuncAnimation(fig1, update_line, 25, 
                                   fargs=(data, data2, l, l2),
                                   interval=50, blit=True)
plt.show()