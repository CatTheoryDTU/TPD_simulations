from ase.thermochemistry import HarmonicThermo
from numpy import array

vibs = array([748.581968,
733.510216,
282.45504,
236.097078,
166.082247
])
vib_energies = vibs / 8065.54429  # convert to eV from cm^-1
#trans_barrier_energy = 0.049313   # eV
#rot_barrier_energy = 0.017675     # eV
#sitedensity = 1.5e15              # cm^-2
#rotationalminima = 6
#symmetrynumber = 1
#mass = 30.07                      # amu
#inertia = 73.149                  # amu Ang^-2

thermo = HarmonicThermo(vib_energies=vib_energies)

F = thermo.get_helmholtz_energy(temperature=298.15)

