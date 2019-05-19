import getopt
import os
import sys


def main():
    g = 10
    dt = 0.01
    m = 1.0
    theta = phi = 0.5
    L = 1.0

    try:
        

def print_usage():
    """Prints usage instructions to stderr and exits."""

    output = "Usage: %s [OPTIONS]\n\n" % os.path.basename(__file__)
    output += "    -h, --help                    "
    output += "prints instructions\n"
    output += "    -v, --verbose                 "
    output += "activates verbose mode\n"
    output += "    -g, --gravity=ACCEL           "
    output += "sets the gravitational acceleration\n"
    output += "    -s, --time-step=STEP          "
    output += "sets the simulation time step\n"
    output += "    -m, --mass=mass      "
    output += "sets the mass of each ball\n"
    output += "    -t, --theta=theta,phi     "
    output += "sets the initial angle of each ball\n"
    output += "    -w, --omega=OMEGA1,OMEGA2     "
    output += "sets the initial angular velocity of each ball\n"
    output += "    -L, --rodlen=LEN      "
    output += "sets the rod length for each ball\n"
    output += "    -verlet     "
    output += "using verlet method\n"
    output += "     -rk4        "
    output += "using fourth order runge-kutta\n"
    sys.stderr.write(output)
    sys.exit(0)

def print_error(errmsg):
    """Prints an error message and exits with an error code (1)."""

    sys.stderr.write("Error: %s\n" % errmsg)
    sys.exit(1)

if __name__ == "__main__":
    main()