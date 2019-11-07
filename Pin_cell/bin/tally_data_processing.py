import matplotlib.pyplot as plt
import numpy as np

import openmc

# Load the statepoint file
sp = openmc.StatePoint('statepoint.20.h5')

tally = sp.get_tally(scores=['flux'])
print(tally)



#standard and mean deviation
print(tally)
(tally.mean, tally.std_dev)

flux = tally.get_slice(scores=['flux'])
#flux =tally.get_slice(scores=['flux'], filters = [openmc.MeshFilter], filter_bins=[-0.62992,0.62992], squeeze = True)
#thermal_flux=flux.get_slice(filters=[openmc.EnergyFilter], filter_bins=[((100., 20.e6),)])
thermal_flux=flux.get_slice(filters=[openmc.EnergyFilter], filter_bins=[((0., 4.),)])

fission = tally.get_slice(scores=['fission'])
thermal_fission=fission.get_slice(filters=[openmc.EnergyFilter], filter_bins=[((0., 4.),)])

nu_fission = tally.get_slice(scores=['nu-fission'])

thermal_flux.std_dev.shape = (100, 100)
thermal_flux.mean.shape = (100, 100)

thermal_fission.std_dev.shape = (100, 100)
thermal_fission.mean.shape = (100, 100)

"""
flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
fission.std_dev.shape = (100, 100)
fission.mean.shape = (100, 100)
"""
"""
fig = plt.subplot(121)
#fig.title = "E neutron flux"
fig.imshow(thermal_flux.mean)
fig2 = plt.subplot(122)
#fig2.title = "Thermal fission events"
fig2.imshow(thermal_fission.mean)
"""
fig = plt.figure()
img = plt.imshow(thermal_flux.mean)
plt.title('Neutron flux, s^(-1)cm^(-2)');
plt.xlabel('x')
plt.ylabel('y')
plt.grid
plt.colorbar(orientation='vertical');
plt.show()

# Determine relative error
relative_error = np.zeros_like(flux.std_dev)
nonzero = flux.mean > 0
relative_error[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]

# distribution of relative errors
ret = plt.hist(relative_error[nonzero], bins=50)
plt.title("Relative error")
plt.show()

#sp.source['E']

# Create log-spaced energy bins from 1 keV to 10 MeV
#energy_bins = np.linspace(0,10000000, 100)
# Create log-spaced energy bins from 1 keV to 10 MeV
energy_bins = np.logspace(-2,7)

# Calculate pdf for source energies
probability, bin_edges = np.histogram(sp.source['E'], energy_bins, density=True)

# Make sure integrating the PDF gives us unity
print(sum(probability*np.diff(energy_bins)))

# Plot source energy PDF
#plt.plot(bin_edges[:-1], probability*np.diff(energy_bins), drawstyle='steps')
plt.semilogx(energy_bins[:-1], probability*np.diff(energy_bins), linestyle='steps')
#plt.title("Neutron energy flux distribution, enrichment 0.1%")
plt.xlabel('Energy (eV)')
plt.ylabel('Probability/eV')
plt.show()



#plt.semilogx(flux.sum, linestyle='steps')
#plt.title("Neutron energy flux distribution, enrichment 0.1%")
#plt.xlabel('Energy (eV)')
#plt.ylabel('Probability/eV')
#plt.show()
#

plt.quiver(sp.source['r']['x'], sp.source['r']['y'],
           sp.source['u']['x'], sp.source['u']['y'],
           np.log(sp.source['E']), cmap='jet', scale=20.0)
plt.colorbar()
plt.xlim((-0.5,0.5))
plt.ylim((-0.5,0.5))


plt.show()


tally_ex = sp.get_tally(name = 'tally 3')
flux_ex = tally_ex.get_slice(scores=['flux'])
print(flux_ex)
print(flux_ex.mean)
#sp.source['flux']

bins = np.logspace(-3,7,49) #Must be the same for the specified tallies
plt.loglog(bins,flux_ex.mean[:,0,0])
plt.xlabel('E, eV')
plt.ylabel('Neutron flux, s^(-1)cm^(-2)')
plt.show()


#plt.quiver(sp.source['r']['x'], sp.source['r']['y'],
#           sp.source['u']['x'], sp.source['u']['y'],
#           np.log(sp.source['E']), cmap='jet', scale=20.0)
#plt.colorbar()
#plt.xlim((-0.5,0.5))
#plt.ylim((-0.5,0.5))


