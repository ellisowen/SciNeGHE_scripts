from gammapy.image import image_profile, lon_lat_rectangle_mask, coordinates
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from gammapy.image.utils import WCS
from astropy.units import Quantity
from gammapy.spectrum.flux_point import _energy_lafferty_power_law
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image import make_empty_image, catalog_image, binary_disk
from gammapy.image.utils import cube_to_image, solid_angle
from gammapy.data import SpectralCube

flux_file = raw_input('Flux Map: ')
exposure_file = raw_input('Exposure Map: ')
spec_ind = input('Spectral Index (for reprojection): ')
flux_hdu = fits.open(flux_file)[1]
flux_wcs = WCS(flux_hdu.header)
energy_flux = Quantity([_energy_lafferty_power_law(10000, 500000, spec_ind)], 'MeV')
flux_data = np.zeros((1, 1800, 3600))
flux_data[0] = Quantity(flux_hdu.data, '')
flux_spec_cube = SpectralCube(data=flux_data, wcs=flux_wcs, energy=energy_flux)

exposure_hdu = fits.open(exposure_file)[0]
exposure_wcs = WCS(exposure_hdu.header)
energy_exp = Quantity([10000, 500000], 'MeV')
exposure_data = Quantity(exposure_hdu.data, '')
exposure_spec_cube = SpectralCube(data=exposure_data, wcs=exposure_wcs, energy=energy_exp)
exposure_spec_cube = exposure_spec_cube.reproject_to(flux_spec_cube)

counts_data = flux_spec_cube.data * exposure_spec_cube.data
counts = fits.ImageHDU(data=counts_data[0], header=flux_hdu.header)



binsz=input('Bin size: ')
label1='Total Flux'
filename1 = flux_file
smooth = 0.1
hdu = fits.open(filename1)[1]
lons, lats = coordinates(hdu)
lat_profile = image_profile(profile_axis='lat', image=hdu, lats=[-10, 10], lons=[-100, 100], binsz=binsz, errors=True, standard_error=0.1, counts=counts)
plt.figure(figsize=(5,4))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
flux = np.float64(lat_profile['BIN_VALUE'])
lats = lat_profile['GLAT_MIN']
profile1 = np.histogram(lats, bins=np.sort(lats), weights=np.float64(flux))
xstep=binsz
y_smooth_1 = gaussian_filter(profile1[0], smooth / xstep)
lats_x = lat_profile['GLAT_MIN'][1:] - 0.5 * binsz
plt.plot(lats_x, y_smooth_1, label='{0}'.format(label1))
plt.hlines(0, -5, 5)
plt.xlabel(r'Galactic Latitude/$deg$', fontsize=10)
plt.ylabel(r'Flux/ph cm$^{-2}$ s$^{-1}$', fontsize=10)
plt.xlim([-5, 5])
plt.tick_params(axis='x', labelsize=10)
plt.grid(b=True, which='major', color='0.75', linewidth=0.5)
plt.legend(prop={'size':8}, loc='lower right')
filename=raw_input('Output filename: ')
plt.savefig(filename)
