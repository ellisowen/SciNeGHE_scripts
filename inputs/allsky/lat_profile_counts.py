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

counts_file = raw_input('Counts Map: ')
background_file = raw_input('Background Map: ')
exposure_file = raw_input('Exposure Map: ')
spec_ind = input('Spectral Index (for reprojection): ')
counts_hdu = fits.open(counts_file)[0]
counts_wcs = WCS(counts_hdu.header)
energy_counts = Quantity([_energy_lafferty_power_law(10000, 500000, spec_ind)], 'MeV')
counts_data = np.zeros((1, 1800, 3600))
counts_data[0] = Quantity(counts_hdu.data, '')
counts_spec_cube = SpectralCube(data=counts_data, wcs=counts_wcs, energy=energy_counts)

exposure_hdu = fits.open(exposure_file)[0]
exposure_wcs = WCS(exposure_hdu.header)
energy_exp = Quantity([10000, 500000], 'MeV')
exposure_data = Quantity(exposure_hdu.data, '')
exposure_spec_cube = SpectralCube(data=exposure_data, wcs=exposure_wcs, energy=energy_exp)
exposure_spec_cube = exposure_spec_cube.reproject_to(counts_spec_cube)

flux_data = counts_spec_cube.data / exposure_spec_cube.data
flux = fits.ImageHDU(data=flux_data[0], header=counts_hdu.header)



binsz=input('Bin size: ')

label1='Total Flux'
smooth = 0.2
lons, lats = coordinates(flux)
lat_profile_total = image_profile(profile_axis='lat', image=flux, lats=[-10, 10], lons=[-100, 100], binsz=binsz, errors=True, counts=counts_hdu)

label2='Separated Background'
background = fits.open(background_file)[1]
lat_profile_background = image_profile(profile_axis='lat', image=background, lats=[-10, 10], lons=[-100, 100], binsz=binsz, errors=True, counts=counts_hdu)

label3='Fermi Diffuse Model'
diffuse = fits.open('diffuse_correct.fits')[1]
lat_profile_diffuse = image_profile(profile_axis='lat', image=diffuse, lats=[-10, 10], lons=[-100, 100], binsz=binsz, errors=False)
plt.clf()
plt.figure(figsize=(7, 6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
flux = np.float64(lat_profile_total['BIN_VALUE'])
background = np.float64(lat_profile_background['BIN_VALUE'])
diffuse_pro = np.float64(lat_profile_diffuse['BIN_VALUE'])
lats = lat_profile_total['GLAT_MIN']
profile1 = np.histogram(lats, bins=np.sort(lats), weights=np.float64(flux))
profile2 = np.histogram(lats, bins=np.sort(lats), weights=np.float64(background))
profile3 = np.histogram(lats, bins=np.sort(lats), weights=np.float64(diffuse_pro))
xstep=binsz
y_smooth_1 = gaussian_filter(profile1[0], smooth / xstep)
y_smooth_2 = gaussian_filter(profile2[0], smooth / xstep)
y_smooth_3 = gaussian_filter(profile3[0], smooth / xstep)
lats_x = lat_profile_total['GLAT_MIN'][1:] - 0.5 * binsz
plt.plot(lats_x, 10 * y_smooth_1, color='grey', label='{0}'.format(label1))
plt.plot(lats_x, 10 * y_smooth_3, color='blue', label='{0}'.format(label3))
plt.plot(lats_x, 10 * y_smooth_2, color='red', label='{0}'.format(label2))
plt.hlines(0, -4.9, 4.9)
plt.xlabel(r'Galactic Latitude/$deg$', fontsize=25)
plt.ylabel(r'Flux/ph cm$^{-2}$ s$^{-1}$ deg$^{-1}$', fontsize=25)
plt.xlim([-4.9, 4.9])
plt.ylim(ymin=0)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.grid(b=True, which='major', color='0.75', linewidth=0.5)
filename=raw_input('Output filename: ')
plt.savefig(filename)
