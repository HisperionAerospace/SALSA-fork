import numpy as np
import scipy as sp
import spiceypy as spice
import matplotlib.pyplot as plt

class options:
    def __init__(self, deltaRA, deltaDec):
        self.deltaRA = deltaRA/180*np.pi     # rad
        self.deltaDec = deltaDec/180*np.pi      # rad

def createInputFilesAtTimeIndex(sat, i):
    orbPosSat = sat.states[0:3, i]
    RASAT, DECSAT = cart2radec(orbPosSat)
    RSAT = np.linalg.norm(orbPosSat)
    SAT_DATA = np.vstack((orbPosSat[0],orbPosSat[1],orbPosSat[2], RASAT, DECSAT, RSAT))
    np.savetxt("INPUT/coord_sat.dat", SAT_DATA, fmt='%.16e', delimiter=',')
    PosSun = sat.sun_vector[0:3, i]
    RASUN, DECSUN = cart2radec(PosSun)
    RSUN = np.linalg.norm(PosSun)
    SUN_DATA = np.vstack((PosSun[0],PosSun[1],PosSun[2], RASUN, DECSUN, RSUN))
    np.savetxt("INPUT/coord_sun.dat", SUN_DATA, fmt='%.16e', delimiter=',')

def cart2radec(vec):
    DEC = np.atan2(vec[2], np.sqrt(vec[0]**2+vec[1]**2))
    RA = np.arctan2(vec[1], vec[0])
    return RA, DEC

def load_kernels():
    # Load all kernels
    spice.kclear()
    spice.furnsh('kernels/mykernel.tm')

def gen_sunAndMoonFile(sat, filename_sun, filename_moon):
    N = sat.times.size
    sun_positions = np.zeros((N, 3))
    moon_positions = np.zeros((N, 3))
    ref = sat.states.transpose()
    for i, t in enumerate(sat.times):
        sun_state = spice.spkezr('SUN', t, 'J2000', 'LT', 'EARTH')
        sun_positions[i, :] = (sun_state[0][0:3])  # Extract position components
        moon_state = spice.spkezr('MOON', t, 'J2000', 'LT', 'EARTH')
        moon_positions[i, :] = (moon_state[0][0:3]-ref[i][0:3])  # Extract position components
    sat.sun_vector = sun_positions.T
    sat.moon_vector = moon_positions.T
    # Stack into 4 columns: time, x, y, z
    table_sun = np.column_stack((sat.times/60.0, sun_positions))
    # Write with space separator, no headers
    np.savetxt(filename_sun, table_sun, fmt='%.16e', delimiter=' ')
    # Stack into 4 columns: time, x, y, z
    table_moon = np.column_stack((sat.times/60.0, moon_positions))
    # Write with space separator, no headers
    np.savetxt(filename_moon, table_moon, fmt='%.16e', delimiter=' ')

def gen_visibleFile(sat, filename, options):
    # Generates the observation map for each time, where all RADEC combinations are observable
    ra = np.linspace(0 + options.deltaRA / 2, 2*np.pi, int(2*np.pi/options.deltaRA))
    dec = np.linspace(-np.pi/2 + options.deltaDec / 2, np.pi/2, int(np.pi/options.deltaDec))
    # Generate the vectors so: each time has all radec combninatios
    ra_repeated = np.repeat(ra, len(dec)) # Repeat each RA for all DECs
    dec_repeated = np.tile(dec, len(ra)) # Tile DECs for all RAs
    sat.ra_targets = ra_repeated
    sat.dec_targets = dec_repeated
    # Generate the table
    table = np.column_stack((ra_repeated, dec_repeated))
    # Write the file
    np.savetxt(filename, table, fmt='%.16e', delimiter=',') 
    n_elem = np.array(len(ra_repeated)).reshape(1)
    np.savetxt("./INPUT/number_of_targets.dat", n_elem, delimiter=',', fmt='%d') 
    return 0

def gen_PST():
    tht = np.linspace(0,np.pi/2,100)
    pst = np.cos(tht)**2
    table = np.column_stack((tht*180/np.pi, pst))
    np.savetxt("DATA/pst.dat", table, fmt='%.16e', delimiter=' ')
    return 0

# Class for the RHS of the TBP
class tbp_rhs:
    
    def __init__(self):
        self.mu = spice.bodvrd('EARTH', 'GM', 1)[1][0]  # Gravitational parameter of Earth in km^3/s^2

    def eval(self, t, x):
        xdot = np.zeros(6)
        xdot[0:3] = x[3:6]
        r = np.linalg.norm(x[0:3])
        xdot[3:6] = -self.mu * x[0:3] / r**3
        return xdot

class Satellite:
    times = np.zeros(1)
    states = np.zeros(6)
    sun_vector = np.zeros(3)
    moon_vector = np.zeros(3)

    def __init__(self, state0):
        self.state0 = state0
        self.rhs = tbp_rhs()

    def get_solution(self, t1, t2, dt=60.0):
        sol = sp.integrate.solve_ivp(
        self.rhs.eval,
        (t1, t2),
        self.state0,
        method="RK45",
        rtol = 1E-12,
        atol = 1E-12,             
        t_eval=np.arange(t1, t2, dt)
        )
        self.times = sol.t
        self.states = sol.y
    def plot_orbit(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")  # 3D axes
        ax.plot(self.states[0], self.states[1], self.states[2], label="3D line")
        ax.legend()
        plt.show()

    def gen_datfile(self, filename):
        # Extract position components (rows 0,1,2)
        x = self.states[0, :]
        y = self.states[1, :]
        z = self.states[2, :]
        
        # Stack into 4 columns: time, x, y, z
        table = np.column_stack((self.times/60.0, x, y, z)) # time in minutes
        
        # Write with space separator, no headers
        np.savetxt(filename, table, fmt='%.16e', delimiter=' ')
        times_rel = self.times/60.0 - self.times[0]/60.0
        orbitid = 1
        table_minutes = np.column_stack((np.full_like(times_rel, orbitid, dtype=np.int32), times_rel, self.times/60.0))
        np.savetxt("DATA/minute_table.dat", table_minutes, fmt='%.16e', delimiter=',')


