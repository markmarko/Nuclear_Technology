import openmc

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 20
inactive = 2
particles = 10000

###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(material_id=1, name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
#uo2.add_element('U', 1., enrichment=1.6)#2.4
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

helium = openmc.Material(material_id=2, name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(material_id=3, name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(material_id=4, name='Borated water (H2O)')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_nuclide('H1', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

heavy_water = openmc.Material(material_id=5, name='Heavy water (D2O)')
heavy_water.set_density('g/cm3', 1.11)
heavy_water.add_nuclide('H2', 2)
heavy_water.add_element('O', 1)
heavy_water.add_s_alpha_beta('c_D_in_D2O')

light_water = openmc.Material(material_id=6, name='Normal water')
light_water.set_density('g/cm3', 0.9998396)
light_water.add_element('H', 2)
light_water.add_element('O', 1)
light_water.add_s_alpha_beta('c_H_in_H2O')

carbon = openmc.Material(material_id=7, name='Carbon')
carbon.set_density('g/cm3', 2.267)
carbon.add_element('C', 1.0)

vapour = openmc.Material(material_id=8, name='Vapour water')
vapour.set_density('g/cm3', 0.0361)
vapour.add_element('H', 2)
vapour.add_element('O', 1)
vapour.add_s_alpha_beta('c_H_in_H2O')

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, helium, zircaloy, borated_water, heavy_water, light_water, carbon, vapour])
materials_file.export_to_xml()

###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate ZCylinder surfaces
factor = 1.0
factor_p = 1.0

entropy_mesh_value = 0.39218			    # 0.39218

fuel_r = 0.39218*factor						# 0.39218
clad_ir = 0.40005*factor					# 0.40005
clar_or = 0.45720*factor					# 0.45720
void_r = 0.5772*factor                     # 0.57720
plane_value = 0.62992*factor_p		        # 0.62992

fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, r=fuel_r, name='Fuel OR')
clad_ir = openmc.ZCylinder(surface_id=2, x0=0, y0=0, r=clad_ir, name='Clad IR')
clad_or = openmc.ZCylinder(surface_id=3, x0=0, y0=0, r=clar_or, name='Clad OR')
void_r= openmc.ZCylinder(surface_id=4, x0=0, y0=0, r=void_r, name='Clad OR')

left = openmc.XPlane(surface_id=5, x0=-plane_value, name='left')
right = openmc.XPlane(surface_id=6, x0=plane_value, name='right')
bottom = openmc.YPlane(surface_id=7, y0=-plane_value, name='bottom')
top = openmc.YPlane(surface_id=8, y0=plane_value, name='top')

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
top.boundary_type = 'reflective'
bottom.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(cell_id=1, name='cell 1')
gap = openmc.Cell(cell_id=2, name='cell 2')
clad = openmc.Cell(cell_id=3, name='cell 3')
void = openmc.Cell(cell_id=4, name='cell 4')
water = openmc.Cell(cell_id=5, name='cell 5')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
gap.region = +fuel_or & -clad_ir
clad.region = +clad_ir & -clad_or
void.region = +clad_or & -void_r
water.region = +void_r & +left & -right & +bottom & -top

# Register Materials with Cells
fuel.fill = uo2
gap.fill = helium
clad.fill = zircaloy

# ********************************
water.fill = light_water
void.fill = light_water
# ********************************

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
root.add_cells([fuel, gap, clad, water, void])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()

###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-plane_value, -plane_value, -1, plane_value, plane_value, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-entropy_mesh_value, -entropy_mesh_value, -1.e50]
entropy_mesh.upper_right = [entropy_mesh_value, entropy_mesh_value, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings_file.entropy_mesh = entropy_mesh
settings_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC tallies.xml file
###############################################################################

tallies_file = openmc.Tallies()

energy_filter = openmc.EnergyFilter([0., 0.625, 20.0e6])

# Instantiate flux Tally in moderator and fuel
tally = openmc.Tally(name='flux')
tally.filters = [openmc.CellFilter([fuel, water])]
tally.filters.append(energy_filter)
tally.scores = ['flux']
tallies_file.append(tally)

# Instantiate reaction rate Tally in fuel
tally = openmc.Tally(name='fuel rxn rates')
tally.filters = [openmc.CellFilter(fuel)]
tally.filters.append(energy_filter)
tally.scores = ['nu-fission', 'scatter']
tally.nuclides = ['U238', 'U235']
tallies_file.append(tally)

# Instantiate reaction rate Tally in moderator
tally = openmc.Tally(name='moderator rxn rates')
tally.filters = [openmc.CellFilter(water)]
tally.filters.append(energy_filter)
tally.scores = ['absorption', 'total']
tally.nuclides = ['O16', 'H1']
tallies_file.append(tally)

# Instantiate a tally mesh
mesh = openmc.RegularMesh(mesh_id=1)
#mesh.dimension = [1, 1, 1]
mesh.dimension = [100, 100, 1]
mesh.lower_left = [-plane_value, -plane_value, -100.]
mesh.width = [2*plane_value, 2*plane_value, 200.]
meshsurface_filter = openmc.MeshSurfaceFilter(mesh)

###############################################################################
#                               Tally arithmetics
###############################################################################
"""
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),openmc.CellFilter([fuel])]
tallies_file.append(fuel_therm_abs_rate)

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
abs_rate = openmc.Tally(name='abs. rate')
fiss_rate.scores = ['nu-fission']
abs_rate.scores = ['absorption']
#tallies_file += (fiss_rate, abs_rate)
tallies_file.append(fiss_rate)
tallies_file.append(abs_rate)

tallies_file.export_to_xml()
"""
###############################################################################
#                                  Plotting
###############################################################################

# Instantiate a tally mesh
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100, 1]
mesh.lower_left = [-plane_value, -plane_value, -1.e50]
mesh.upper_right = [plane_value, plane_value, 1.e50]

# Instantiate some tally Filters
energy_filter = openmc.EnergyFilter([0., 4., 100., 20.e6])
mesh_filter = openmc.MeshFilter(mesh)

# Instantiate the Tally
tally = openmc.Tally(tally_id=1, name='tally 1')
tally.filters = [energy_filter, mesh_filter]
tally.scores = ['flux', 'fission', 'nu-fission']

# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies([tally])
#tallies_file.append(tally)
tallies_file.export_to_xml()
 


plot1 = openmc.Plot()
plot1.filename = 'materials-xy'
plot1.origin = [0, 0, 0]
plot1.basis = 'xy'
plot1.width = [1.0, 1.0]
plot1.pixels = [3000, 3000]
plot1.color = "material"

# Plot inline in jupyter or QtConsole
#openmc.plot_inline(plot1)
plots = openmc.Plots([plot1])
plots.export_to_xml()
#openmc.plot_geometry()

openmc.run()

"""
###############################################################################
#                 			Post-Processing
###############################################################################

# Load the statepoint file
sp = openmc.StatePoint('statepoint.10.h5')

tally = sp.get_tally(scores=['flux'])
print(tally)

#tally.sum

#standard and mean deviation
print(tally.mean.shape)
(tally.mean, tally.std_dev)

flux = tally.get_slice(scores=['flux'])
fission = tally.get_slice(scores=['fission'])
print(flux)

flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
fission.std_dev.shape = (100, 100)
fission.mean.shape = (100, 100)

fig = plt.subplot(121)
fig.imshow(flux.mean)
fig2 = plt.subplot(122)
fig2.imshow(fission.mean)

# Determine relative error
relative_error = np.zeros_like(flux.std_dev)
nonzero = flux.mean > 0
relative_error[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]

# distribution of relative errors
ret = plt.hist(relative_error[nonzero], bins=50)

#sp.source['E']

# Create log-spaced energy bins from 1 keV to 10 MeV
energy_bins = np.logspace(3,7)

# Calculate pdf for source energies
probability, bin_edges = np.histogram(sp.source['E'], energy_bins, density=True)

# Make sure integrating the PDF gives us unity
print(sum(probability*np.diff(energy_bins)))

# Plot source energy PDF
plt.semilogx(energy_bins[:-1], probability*np.diff(energy_bins), linestyle='steps')
plt.xlabel('Energy (eV)')
plt.ylabel('Probability/eV')

plt.quiver(sp.source['r']['x'], sp.source['r']['y'],
           sp.source['u']['x'], sp.source['u']['y'],
           np.log(sp.source['E']), cmap='jet', scale=20.0)
plt.colorbar()
plt.xlim((-0.5,0.5))
plt.ylim((-0.5,0.5))

"""

