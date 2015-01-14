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

    print len(optimized)
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