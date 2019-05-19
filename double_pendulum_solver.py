import math, numpy

class DoublePendulum:
    def __init__(self, g, m, theta, phi, L, w_1, w_2):
        """
        g - The gravitational acceleration
        m - the mass of two balls
        theta - initial angle of ball #1
        phi - intial angle of ball #2
        L - length of the rod for both balls
        """

        self.g = g
        self.m = m
        self.theta = theta 
        self.phi = phi
        self.L = L

        self.w_1 = w_1 # initial angular velocity ball 1
        self.w_2 = w_2 # intial angular velocity ball 2

        self.nat_w = self.g/self.L

        self.xi = [self.theta, self.w_1, self.phi, self.w_2] #y(0)
    
    def compute_right_hand_size(self, theta, w_1, phi, w_2):
        """
        theta - the angle of first ball
        phi - the angle of second ball
        w_1 - angular velocity of first ball
        w_2 - angular velocity of second ball
        """

        g_numerator_1 = 2 * (self.nat_w)**2 ) * ( math.sin(theta) * math.cos(theta - phi) - math.sin(phi) )
        g_numerator_2 = ((math.cos(theta - phi) * (w_2**2)) + 2 * (w_1**2)) * math.sin(theta - phi)
        g_numerator = g_numerator_1 + g_numerator_2
        g_denominator = 1 + (math.sin(theta -phi))**2
        g = g_numerator / g_denominator 

        f_numerator_1 = (self.nat_w ** 2) * ( math.sin(phi) * math.cos(theta - phi) - 2 * math.phi(theta))
        f_numerator_2 = (w_2 ** 2 + w_1 ** 2 * math.cos(theta - phi)) * math.sin(theta - phi)
        f_numerator = f_numerator_1 - f_numerator_2
        f_denominator = 1 + (math.sin(theta - phi)) ** 2
        f = f_numerator / f_denominator

        return numpy.array([w_1, f, w_2, g])

    def integrate_verlet(self, dt):
        #advance one time step
        #verlet integration
        #temp copies of current model's variables 
        theta = self.theta
        phi = self.phi
        w_1 = self.w_1
        w_2 = self.w_2

        #calculate y(k*dt)
        self.xi = self.xi + self.compute_right_hand_size(theta, w_1, phi, w_2) * dt

        self.theta = self.xi[0]
        self.w_1 = self.xi[1]
        self.phi = self.xi[2]
        self.w_2 = self.xi[3]

    def integrate_runge_kutta(self, dt):
        #advance one time step
        #runge kutta integration

        inputs = numpy.array([self.theta, self.w_1, self.phi, self.w_2])

        k_1 = self.compute_right_hand_size(*inputs)
        k_2 = self.compute_right_hand_size(*(inputs + dt * k_1 / 2))
        k_3 = self.compute_right_hand_size(*(inputs + dt * k_2 / 2))
        k_4 = self.compute_right_hand_size(*(inputs + dt * k_3 ))

        self.xi = self.xi + (k_1 + 2*k_2 + 2*k_3 + k_4) * dt / 6

        self.theta = self.xi[0]
        self.w_1 = self.xi[1]
        self.phi = self.xi[2]
        self.w_2 = self.xi[3]

        



    