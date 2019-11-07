import matplotlib.pyplot as plt
import numpy as np
import openmc
from os import listdir
from os.path import isfile, join
import openmc.deplete

#Creating a particle filter
sp = openmc.StatePoint('statepoint.20.h5')

sp_tally = sp.get_tally(scores=['flux','fission','nu-fission'])

flux = sp_tally.get_slice(scores=['flux'])
neutron_flux = flux.get_slice(filters = [openmc.ParticleFilter], filter_bins=[(('neutron'),)], squeeze = True)

fission = sp_tally.get_slice(scores=['fission'])
neutron_fission = fission.get_slice(filters = [openmc.ParticleFilter],filter_bins = [(('neutron'),)], squeeze = True)



neutron_flux.std_dev.shape=(200,200)
neutron_flux.mean.shape=(200,200)

neutron_fission.std_dev.shape=(200,200)
neutron_fission.mean.shape=(200,200)



fig = plt.subplot(121)
fig.imshow(neutron_flux.mean)
fig2 = plt.subplot(122)
fig2.imshow(neutron_fission.mean)
plt.show()

tally_ex = sp.get_tally(name = 'tally 4')
flux_ex = tally_ex.get_slice(scores=['flux'])
print(flux_ex)
print(flux_ex.mean)
#sp.source['flux']

bins = np.logspace(-3,7,49) #Must be the same for the specified tallies
plt.loglog(bins,flux_ex.mean[:,0,0])
plt.xlabel('E, eV')
plt.ylabel('Neutron flux, s^(-1)cm^(-2)')
plt.show()
