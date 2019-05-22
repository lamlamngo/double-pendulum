import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import double_pendulum_solver as dps

def main():
    #default values
    G = 9.8
    L = 1.0
    M = 1.0
    dt = 0.025

    #initial conditions
    theta = np.radians(10.0)
    phi = np.radians(10.0)
    w_1 = np.radians(0.0)
    w_2 = np.radians(0.0)

    inputs = np.array([G, M, theta, phi, L , w_1, w_2, dt])
    dp = dps.DoublePendulum(*inputs)

    thetas = []
    thetas.append(theta)
    w_1s = []
    w_1s.append(w_1)

    fig, (ax, ax1) = plt.subplots(1,2)
    plt.subplots_adjust(wspace = 0.36)
    ax.set(xlabel='x (m)', ylabel='y (m)')
    ax1.set(xlabel='θ (rad)', ylabel='ξ (rad/s)')

    ax.grid()
    ax1.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    line1, = ax1.plot(thetas, w_1s, '-', lw=0.5)

    time_template = 'time = %.2fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    ax1.set_xlim(-abs(theta) - 0.1, abs(theta) + 0.1)
    ax1.set_ylim(-abs(w_1) - 0.1, abs(w_1) + 0.1)

    ax.set_xlim(-2*L - 0.1, 2*L + 0.1)
    ax.set_ylim(-2*L - 0.1, 2*L + 0.1)

    def init():
        line.set_data([],[])
        line1.set_data(thetas, w_1s)
        time_text.set_text('')
        return line, line1, time_text
    
    def animate(i):
        dp.integrate_runge_kutta(dt)

        current_max_x = ax1.get_xlim()[1]
        current_max_y = ax1.get_ylim()[1]

        current_max_x = abs(dp.theta) + 0.1 if abs(dp.theta) > current_max_x else current_max_x
        current_max_y = abs(dp.w_1) + 0.1 if abs(dp.w_1) > current_max_y else current_max_y

        ax1.set_xlim(-current_max_x, current_max_x)
        ax1.set_ylim(-current_max_y, current_max_y)

        x1 = L * math.sin(dp.theta)
        x2 = x1 + L * math.sin(dp.phi)
        y1 = -L * math.cos(dp.theta)
        y2 = y1 - L * math.cos(dp.phi)

        thetas.append(dp.theta)
        w_1s.append(dp.w_1)

        line1.set_data(thetas, w_1s)
        line.set_data([0,x1,x2], [0,y1,y2])
        time_text.set_text(time_template % (i*dt))

        return line, line1, time_text
    
    ani = animation.FuncAnimation(fig, animate, interval= dt*1000, blit=False, init_func=init)
    plt.show()

if __name__ == "__main__":
    main()