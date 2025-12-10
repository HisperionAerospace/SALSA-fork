import numpy as np
import inbuilt.orbit_comp as orb

SolutionMatrix = np.load("OUTPUT/SolutionMatrix.npy")
IndexMap = np.load("OUTPUT/IndexMap.npy", allow_pickle=True).item()
sat = np.load("OUTPUT/sat.npy", allow_pickle=True).item()
# To access the solution at a pair RA,DEC at time index t: SolutionMatrix[t, IndexMap[(RA, DEC)]]
RA_VEC = sat.ra_targets[0]  # Example RA in radians
DEC_VEC = sat.dec_targets[0]  # Example DEC in radians
plot_points = SolutionMatrix[:, IndexMap[(RA_VEC, DEC_VEC)]]
i = 0
for point in plot_points:
    print("%3.1f %.16f" % (i,point))
    i=i+1

import matplotlib.pyplot as plt
plt.plot(sat.times/60.0, plot_points)
plt.xlabel("Time [min]")
plt.ylabel("Stray Light Level")
plt.title("Stray Light Level at RA=%.2f deg, DEC=%.2f deg" % (RA_VEC*180/np.pi, DEC_VEC*180/np.pi))
plt.grid()
plt.show()