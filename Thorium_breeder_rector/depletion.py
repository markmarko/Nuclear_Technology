import matplotlib.pyplot as plt
import numpy as np
import openmc
from os import listdir
from os.path import isfile, join
import openmc.deplete

chain = openmc.deplete.Chain.from_xml("./chain_casl.xml")
print(chain.nuclide_dict)

results = openmc.deplete.ResultsList.from_hdf5("./depletion_results.h5")

time, k = results.get_eigenvalue()

time /= (24 * 60 * 60)

plt.errorbar(time, k[:, 0], yerr=k[:, 1])
plt.xlabel("Time [d]")
plt.ylabel("$k_{eff}\pm \sigma$");

plt.show()

time, u3 = results.get_atoms("11", "U233")
time, xe135 = results.get_atoms("11", "Xe135")
time, th232 = results.get_atoms("11", "Th232")
time, pu239 = results.get_atoms("11", "Pu239")
time, pa231 = results.get_atoms("11", "Pa231")
time, u232 = results.get_atoms("11", "U232")

time /= (24 * 60 * 60)

plt.plot(time, u3, label="U233")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - U233");

plt.show()

plt.plot(time, xe135, label="Xe135")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - Xe135");

plt.show()

plt.plot(time, th232, label="Th232")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - Th232");

plt.show()

plt.plot(time, pu239, label="Pu239")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - Pu239");

plt.show()

plt.plot(time, pa231, label="Pa231")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - Pa231");

plt.show()

plt.plot(time, u232, label="U232")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - U232");

plt.show()