import openmc
uo2 = openmc.Material(1, "uo2")
print(uo2)
mat = openmc.Material()
print(mat)

#help(uo2.add_nuclide)
# Add nuclides to uo2
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238', 0.97)
uo2.add_nuclide('O16', 2.0)

uo2.set_density('g/cm3', 10.0)

zirconium = openmc.Material(2, "zirconium")
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.6)

water = openmc.Material(3, "h2o")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)

water.add_s_alpha_beta('c_H_in_H2O')

mats = openmc.Materials([uo2, zirconium, water])

mats = openmc.Materials()
mats.append(uo2)
mats += [zirconium, water]
isinstance(mats, list)

mats.export_to_xml()

water.remove_nuclide('O16')
water.add_element('O', 1.0)

mats.export_to_xml()

print('    ...')

uo2_three = openmc.Material()
uo2_three.add_element('U', 1.0, enrichment=3.0)
uo2_three.add_element('O', 2.0)
uo2_three.set_density('g/cc', 10.0)

sph = openmc.Sphere(r=1.0)

inside_sphere = -sph
outside_sphere = +sph

print((0,0,0) in inside_sphere, (0,0,2) in inside_sphere)
print((0,0,0) in outside_sphere, (0,0,2) in outside_sphere)

z_plane = openmc.ZPlane(z0=0)
northern_hemisphere = -sph & +z_plane

northern_hemisphere.bounding_box

cell = openmc.Cell()
cell.region = northern_hemisphere

# or...
cell = openmc.Cell(region=northern_hemisphere)
cell.fill = water


universe = openmc.Universe()
universe.add_cell(cell)

# this also works
universe = openmc.Universe(cells=[cell])

universe.plot(width=(2.0, 2.0))



