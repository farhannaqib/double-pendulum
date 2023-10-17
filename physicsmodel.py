import matplotlib.pyplot as plt
import numpy as np
import math

import matplotlib.animation as animation

fig, ax = plt.subplots()
g = -9.81
r = 5

theta = 0.5
thetadot = 0

deltat = 0.1
tmax = 10

ax.set(xlim = (-r*1.1, r*1.1), ylim=(-r*1.1, .25))
graph = ax.plot([r*math.sin(theta), 0], [-r*math.cos(theta), 0])[0]

def update(frame):
    global theta, thetadot

    graph.set_xdata([r*math.sin(theta), 0])
    graph.set_ydata([-r*math.cos(theta), 0])
    theta += thetadot*deltat
    thetadot += (g*math.sin(theta)/r)*deltat
    return graph

ani = animation.FuncAnimation(fig=fig, func=update, frames= (int) (tmax/deltat), interval=deltat)
plt.show()