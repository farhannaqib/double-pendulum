import matplotlib.pyplot as plt
import numpy as np
import math

import matplotlib.animation as animation

from dataclasses import dataclass

class DoublePendulum:

    g = 9.81

    @dataclass
    class AttachedMass:
        mass: float
        length: float
        theta: float
        thetadot: float

    def __init__(self, m1, l1, theta1, thetadot1, m2, l2, theta2, thetadot2):
        self.p1 = self.AttachedMass(m1, l1, theta1, thetadot1)
        self.p2 = self.AttachedMass(m2, l2, theta2, thetadot2)

    def A(self):
        return (self.p2.mass * self.p2.length)/((self.p2.mass + self.p1.mass) * self.p1.length) * math.cos(self.p1.theta - self.p2.theta)
    
    def B(self):
        return self.p1.length/self.p2.length * math.cos(self.p1.theta - self.p2.theta)
    
    def C(self): 
        return -(self.p2.mass * self.p2.length) / ((self.p1.mass + self.p2.mass) * self.p1.length) * (self.p2.thetadot ** 2) * math.sin(self.p1.theta - self.p2.theta) - self.g * math.sin(self.p1.theta) / self.p1.length
    
    def D(self):
        return (self.p1.length * (self.p1.thetadot ** 2) * math.sin(self.p1.theta - self.p2.theta) - self.g * math.sin(self.p2.theta)) / self.p2.length
    
    def det(self):
        return 1 - self.A() * self.B()
    
    def updatePoints(self, deltaT): 
        self.p1.theta += self.p1.thetadot * deltaT
        self.p2.theta += self.p2.thetadot * deltaT
        self.p1.thetadot += (self.C() - self.A() * self.D()) / self.det() * deltaT
        self.p2.thetadot += (self.D() - self.B() * self.C()) / self.det() * deltaT

    def getXs(self):
        return [0, self.p1.length * math.sin(self.p1.theta), self.p1.length * math.sin(self.p1.theta) + self.p2.length * math.sin(self.p2.theta)]
    
    def getYs(self):
        return [0, -self.p1.length * math.cos(self.p1.theta), -self.p1.length * math.cos(self.p1.theta) - self.p2.length * math.cos(self.p2.theta)]


fig, ax = plt.subplots()

theta1 = 1
theta1dot = 0
theta2 = -1
theta2dot = 0

l1 = 1
l2 = 1

m1 = 1
m2 = 1

deltat = 0.05
tmax = 10

pendulum = DoublePendulum(m1, l1, theta1, theta1dot, m2, l2, theta2, theta2dot)
ax.set(xlim = (-(l1+l2)*1.1, (l1+l2)*1.1), ylim=(-(l1+l2)*1.1, .25))
graph = ax.plot(pendulum.getXs(), pendulum.getYs())[0]

def update(frame):
    graph.set_xdata(pendulum.getXs())
    graph.set_ydata(pendulum.getYs())
    pendulum.updatePoints(deltat)
    return graph

ani = animation.FuncAnimation(fig=fig, func=update, frames= (int) (tmax/deltat), interval=deltat)
plt.show()