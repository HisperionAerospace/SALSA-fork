import inbuilt.orbit_comp as orb
import inbuilt.compiler as comp
import numpy as np
import os
import sys
import subprocess

# Settings
initial_state = np.array([8000.0, 0.0, 0.0, 0.0, np.sqrt(398600.4355070226/8000), 0.0])  # Initial state vector [km, km/s]
t0 = 0.0  # Initial time in seconds
tf = 24*3600  # Final time in seconds
options = orb.options(360,180)  # Example options for visibility calculations

orb.load_kernels()
sat = orb.Satellite(initial_state)
sat.get_solution(t0,tf, 60) # get solution with time step of 60 seconds
# sat.plot_orbit()
# Generate all necessary input files
sat.gen_datfile("DATA/ORBIT/orbit.dat")
orb.gen_sunAndMoonFile(sat, "DATA/ORBIT/sun.dat", "DATA/ORBIT/moon.dat")
orb.gen_visibleFile(sat, "INPUT/coord_targets.dat", options)
orb.gen_PST()

# compile the code
comp.compiler_core()

# Loop for all the time instants
path = os.getcwd()
SolutionMatrix = np.zeros((len(sat.times),len(sat.ra_targets)))
IndexMap = {(float(x), float(y)): i for i, (x, y) in enumerate(zip(sat.ra_targets, sat.dec_targets))}
for i in range(len(sat.times)):
    # Create the input files for each time
    orb.createInputFilesAtTimeIndex(sat,i)
    os.chdir(os.path.join(os.path.abspath(sys.path[0]), '%s/FORTRAN_CORE_CODE/' % path))
    subprocess.call(["./stray_light"])
    os.chdir(path)
    data = np.loadtxt("./OUTPUT/straylight.out").reshape(-1,3)
    os.remove("./OUTPUT/straylight.out")
    if data is not None and data.size > 0:
        indices = np.array([IndexMap[(x, y)] for x, y in zip(data[:,0], data[:,1])])
        SolutionMatrix[i, indices] = data[:,2]
    if (i) % 10 == 0 or i == len(sat.times)-1:
        print("Completed time index %d / %d" % (i+1, len(sat.times)))

np.save("OUTPUT/SolutionMatrix.npy", SolutionMatrix)
np.save("OUTPUT/IndexMap.npy", IndexMap)
np.save("OUTPUT/sat.npy", sat)




