import getopt
import os
import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
import double_pendulum_solver as dp
import math

def main():
    #default values
    G = 9.8
    L = 1.0
    M = 1.0
    dt = 0.025

    #time stepping range 0 -> 20 seconds, by dt. retyrns list
    time_steps = np.arange(0.0, 20, dt)

    #default initial
    theta = np.radians(120.0)
    phi = np.radians(-10.0)
    w_1 = np.radians(0.0)
    w_2 = np.radians(0.0)

    inputs = np.array([G, M, theta, phi, L, w_1, w_2, dt])
    double_pendulum_system = dp.DoublePendulum(*inputs)

    thetas = []
    phis = []

    thetas.append(theta)
    phis.append(phi)

    i = 1
    #time stepping loop to create values
    while (i < len(time_steps)):
        double_pendulum_system.integrate_runge_kutta(dt)
        thetas.append(double_pendulum_system.theta)
        phis.append(double_pendulum_system.phi)
        i += 1
    
    x_1 = []
    y_1 = []
    x_2 = [] 
    y_2 = []

    #calculate coords
    x_1 = [L * math.sin(theta) for theta in thetas]
    y_1 = [-L * math.cos(theta) for theta in thetas]
    x_2 = [L*math.sin(phis[i]) + x_1[i] for i in range(len(thetas))]
    y_2 = [-L*math.cos(phis[i]) + y_1[i] for i in range(len(thetas))]

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2,2))
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.2fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(i):
        thisx = [0, x_1[i], x_2[i]]
        thisy = [0, y_1[i], y_2[i]]

        line.set_data(thisx, thisy)
        time_text.set_text(time_template % (i*dt))
        return line, time_text   

    interval = 1000 * dt
    
    ani = animation.FuncAnimation(fig, animate, range(0, len(thetas)), interval= interval, blit=True, init_func=init)
    plt.show()

if __name__ == "__main__":
    main()