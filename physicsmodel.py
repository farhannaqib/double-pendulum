import matplotlib.pyplot as plt
import numpy as np
import math

import matplotlib.animation as animation

from dataclasses import dataclass

class DoublePendulum:

    G = 9.81

    @dataclass
    class AttachedMass:
        mass: float
        length: float
        theta: float
        thetadot: float

    def __init__(self, m1, l1, theta1, thetadot1, m2, l2, theta2, thetadot2):
        self.p1 = self.AttachedMass(m1, l1, theta1, thetadot1)
        self.p2 = self.AttachedMass(m2, l2, theta2, thetadot2)
        
    def updatePoints(self, deltaT): 
        # Approximation is more unstable when using old values for theta
        self.p1.theta += self.p1.thetadot * deltaT
        self.p2.theta += self.p2.thetadot * deltaT

        m1 = self.p1.mass
        m2 = self.p2.mass
        l1 = self.p1.length
        l2 = self.p2.length
        theta1 = self.p1.theta
        theta2 = self.p2.theta
        thetadot1 = self.p1.thetadot
        thetadot2 = self.p2.thetadot

        a = (m2 * l2) / ((m1 + m2) * l1) * math.cos(theta1 - theta2)
        b = l1 / l2 * math.cos(theta1 - theta2)
        c = -(m2 * l2) / ((m1 + m2) * l1) * (thetadot2 ** 2) * math.sin(theta1 - theta2) - self.G * math.sin(theta1) / l1
        d = (l1 * (thetadot1 ** 2) * math.sin(theta1 - theta2) - self.G * math.sin(theta2)) / l2

        self.p1.thetadot += (c - a * d) / (1 - a * b) * deltaT
        self.p2.thetadot += (d - b * c) / (1 - a * b) * deltaT

    def getCoord(self, massPoint = AttachedMass):
        return [massPoint.length * math.sin(massPoint.theta), -massPoint.length * math.cos(massPoint.theta)]

    def getXs(self):
        return [0, self.getCoord(self.p1)[0], self.getCoord(self.p1)[0] + self.getCoord(self.p2)[0]]
    
    def getYs(self):
        return [0, self.getCoord(self.p1)[1], self.getCoord(self.p1)[1] + self.getCoord(self.p2)[1]]
    
    def K(self):
        m1 = self.p1.mass
        m2 = self.p2.mass
        l1 = self.p1.length
        l2 = self.p2.length
        theta1 = self.p1.theta
        theta2 = self.p2.theta
        thetadot1 = self.p1.thetadot
        thetadot2 = self.p2.thetadot
        return 1/2.0 * m1 * (l1 * thetadot1) ** 2 + 1/2.0 * m2 * ((l1 * thetadot1) ** 2 + 
                (l2 * thetadot2) ** 2 + (2 * l1 * l2 * thetadot1 * thetadot2 * math.cos(theta1 - theta2)))

    def U(self):
        m1 = self.p1.mass
        m2 = self.p2.mass
        l1 = self.p1.length
        l2 = self.p2.length
        theta1 = self.p1.theta
        theta2 = self.p2.theta
        return -self.G * ((m1 + m2) * l1 * math.cos(theta1) + (m2 * l2 * math.cos(theta2)))

fig, ax = plt.subplots()

theta1 = 0.2
thetadot1 = 0
theta2 = 0.4
thetadot2 = 0

l1 = 1
l2 = 1

m1 = 1
m2 = 1

deltat = 0.05
tmax = 10

pendulum = DoublePendulum(m1, l1, theta1, thetadot1, m2, l2, theta2, thetadot2)
ax.set(xlim = (-(l1+l2)*1.1, (l1+l2)*1.1), ylim=(-(l1+l2)*1.1, .25))
graph = ax.plot(pendulum.getXs(), pendulum.getYs())[0]

def update(frame):
    graph.set_xdata(pendulum.getXs())
    graph.set_ydata(pendulum.getYs())
    pendulum.updatePoints(deltat)
    return graph

ani = animation.FuncAnimation(fig=fig, func=update, frames= (int) (tmax/deltat), interval=deltat)
plt.show()