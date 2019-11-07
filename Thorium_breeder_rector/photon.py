import matplotlib.pyplot as plt
import numpy as np
import openmc
from os import listdir
from os.path import isfile, join
import openmc.deplete

sp = openmc.StatePoint('statepoint.10.h5')

tally = sp.get_tally(scores=['flux'],name = 'flux_xy')
tally.sum
(tally.mean, tally.std_dev)
flux_xy = tally.get_slice(scores=['flux'],filters = [openmc.ParticleFilter], filter_bins = [('neutron','photon')])
neutron_flux_xy = flux_xy.get_slice(filters = [openmc.ParticleFilter], filter_bins = [('neutron',)], squeeze=True)
photon_flux_xy = flux_xy.get_slice(filters = [openmc.ParticleFilter], filter_bins = [('photon',)], squeeze=True)

tally = sp.get_tally(scores=['flux'],name = 'flux_xz')

tally.sum
#print(tally.mean.shape)
(tally.mean, tally.std_dev)
flux_xz = tally.get_slice(scores=['flux'],filters = [openmc.ParticleFilter], filter_bins = [('neutron','photon')])
#fission = tally.get_slice(scores=['fission'])

neutron_flux_xz = flux_xz.get_slice(filters = [openmc.ParticleFilter], filter_bins = [('neutron',)], squeeze=True)
photon_flux_xz = flux_xz.get_slice(filters = [openmc.ParticleFilter], filter_bins = [('photon',)], squeeze=True)

print(flux_xy)
print(photon_flux_xy)
print(neutron_flux_xy)

print(flux_xz)
print(photon_flux_xz)
print(neutron_flux_xz)

neutron_flux_xy.std_dev.shape = (200, 200)
neutron_flux_xy.mean.shape = (200, 200)


photon_flux_xy.std_dev.shape = (200, 200)
photon_flux_xy.mean.shape = (200, 200)

fig = plt.figure()
plt.figure(figsize=(15,15));
img = plt.imshow(neutron_flux_xy.mean,interpolation='spline16');
plt.title('Neutron flux');
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.colorbar( orientation='vertical')
plt.show()

fig = plt.figure()
plt.figure(figsize=(15,15))
img = plt.imshow(photon_flux_xy.mean,interpolation='spline16')
plt.title('Photon flux');
plt.xlabel('x')
plt.ylabel('y')
plt.grid
plt.colorbar(orientation='vertical');
plt.show()