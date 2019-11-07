# Load the statepoint file
import openmc

sp = openmc.StatePoint('statepoint.20.h5')

# Get the fission and absorption rate tallies
fiss_rate = sp.get_tally(name='fiss. rate')
abs_rate = sp.get_tally(name='abs. rate')

keff = fiss_rate / (abs_rate)
print(keff.mean,"+-",keff.std_dev)

# Compute fast fission factor factor using tally arithmetic
therm_abs_rate = sp.get_tally(name='therm. abs. rate')
fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
therm_util = fuel_therm_abs_rate / therm_abs_rate

print(therm_util.mean,"+-",therm_util.std_dev)