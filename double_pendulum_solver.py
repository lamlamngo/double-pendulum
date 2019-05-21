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
        #comment: We actually don't really need the masses of the pendula, as long as they are equal. An interesting feature of the simple 
        #pendulum is that the angular frequency does not depend on the mass!

        self.g = g
        self.m = m
        self.theta = theta 
        self.phi = phi
        self.L = L

        self.w_1 = w_1 # initial angular velocity ball 1
        self.w_2 = w_2 # intial angular velocity ball 2
        
        self.nat_w = self.g/self.L
        
        #define lists for Verlet integration
        self.psi = numpy.array([self.theta, self.phi])
        self.omega = numpy.array([self.w_1, self.w_2])
        
        #Initial value for Psi((n-1)*dt)
        theta_minus1 = self.theta - w_1*dt
        phi_minus1 = self.theta - w_2*dt
        self.psi_minus1 = numpy.([theta_minus1, phi_minus1])

    def compute_right_hand_size(self, theta, w_1, phi, w_2):
        """
        theta - the angle of first ball
        phi - the angle of second ball
        w_1 - angular velocity of first ball
        w_2 - angular velocity of second ball
        """

        g_numerator_1 = 2 * ((self.nat_w) ) * ( math.sin(theta) * math.cos(theta - phi) - math.sin(phi) )
        g_numerator_2 = ((math.cos(theta - phi) * (w_2**2)) + 2 * (w_1**2)) * math.sin(theta - phi)
        g_numerator = g_numerator_1 + g_numerator_2
        g_denominator = 1 + (math.sin(theta -phi))**2
        g = g_numerator / g_denominator 

        f_numerator_1 = (self.nat_w) * ( math.sin(phi) * math.cos(theta - phi) - 2 * math.sin(theta))
        f_numerator_2 = ((w_2 ** 2) + (w_1 ** 2) * math.cos(theta - phi)) * math.sin(theta - phi)
        f_numerator = f_numerator_1 - f_numerator_2
        f_denominator = 1 + (math.sin(theta - phi)) ** 2
        f = f_numerator / f_denominator

        return numpy.array([w_1, f, w_2, g])
   
  def accelerations(self, theta, w_1, phi, w_2):
        """Verlet is a second-order algorithm, which requires the second-order term in the Taylor expansion, namely the acceleration. 
        Therefore, we better work with equations 2.2.7 and 2.2.8 in my note than equation 2.2.12"""
        
        g_numerator_1 = 2 * ((self.nat_w) ) * ( math.sin(theta) * math.cos(theta - phi) - math.sin(phi) )
        g_numerator_2 = ((math.cos(theta - phi) * (w_2**2)) + 2 * (w_1**2)) * math.sin(theta - phi)
        g_numerator = g_numerator_1 + g_numerator_2
        g_denominator = 1 + (math.sin(theta -phi))**2
        g = g_numerator / g_denominator 

        f_numerator_1 = (self.nat_w) * ( math.sin(phi) * math.cos(theta - phi) - 2 * math.sin(theta))
        f_numerator_2 = ((w_2 ** 2) + (w_1 ** 2) * math.cos(theta - phi)) * math.sin(theta - phi)
        f_numerator = f_numerator_1 - f_numerator_2
        f_denominator = 1 + (math.sin(theta - phi)) ** 2
        f = f_numerator / f_denominator
        
        return numpy.array([f,g]) #This is practically identical to the previous method, only now we need a two-entry array

    def integrate_verlet(self, dt):
        #advance one time step
        #verlet integration
        #temp copies of current model's variables 
        theta = self.theta
        phi = self.phi
        w_1 = self.w_1
        w_2 = self.w_2
        clone_Psi = self.psi
  
        #calculate y(k*dt)
        self.psi = self.psi*2 - self.psi_minus1 + self.accelerations(theta, w_1, phi, w_2) * dt 
        self.omega = (self.Psi - clone_Psi)/dt

        self.theta = self.psi[0]
        self.w_1 = self.omega[0]
        self.phi = self.osi[1]
        self.w_2 = self.omega[1]
        self.psi_minus1 = clone_Psi
        
    def integrate_runge_kutta(self, dt):
        #advance one time step
        #runge kutta integration

        inputs = numpy.array([self.theta, self.w_1, self.phi, self.w_2])

        k_1 = self.compute_right_hand_size(*inputs)
        k_2 = self.compute_right_hand_size(*(inputs + dt * k_1 * 0.5))
        k_3 = self.compute_right_hand_size(*(inputs + dt * k_2 * 0.5))
        k_4 = self.compute_right_hand_size(*(inputs + dt * k_3 ))

        delta = (k_1 + 2*k_2 + 2*k_3 + k_4) * dt * 1.0 / 6.0

        self.theta += delta[0]
        self.w_1 += delta[1]
        self.phi += delta[2]
        self.w_2 += delta[3]

        



    
