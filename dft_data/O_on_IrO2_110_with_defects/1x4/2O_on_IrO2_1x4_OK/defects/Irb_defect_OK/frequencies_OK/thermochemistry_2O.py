from ase.thermochemistry import HarmonicThermo
from numpy import array

vibs = array([774.855024,
762.888486,
261.446313,
226.838034,
183.716502
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

