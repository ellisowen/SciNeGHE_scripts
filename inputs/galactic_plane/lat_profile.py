from gammapy.image import image_profile
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

binsz=input('Bin size: ')
label1=raw_input('Type 1: ')
filename1 = raw_input('Filename 1: ')
smooth = 0.2
hdu = fits.open(filename1)[1]
lat_profile = image_profile(profile_axis='lat', image=hdu, lats=[-10, 10], lons=[-100, 100], binsz=binsz, errors=True, standard_error=0.1)
plt.figure(figsize=(5,4))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
flux = np.float64(lat_profile['BIN_VALUE'])
lats = lat_profile['GLAT_MIN']#np.float64(0.5*(lat_profile['GLAT_MIN']+lat_profile['GLAT_MAX']))
profile1 = np.histogram(lats, bins=np.sort(lats), weights=np.float64(flux))
xstep=binsz
y_smooth_1 = gaussian_filter(profile1[0], smooth / xstep)
lats_x = lat_profile['GLAT_MIN'][1:] - 0.5 * binsz#np.float64(0.5*(lat_profile['GLAT_MIN'][1:]+lat_profile['GLAT_MAX'][1:]))
plt.plot(lats_x, y_smooth_1, label='{0}'.format(label1))
plt.hlines(0, -5, 5)
plt.xlabel(r'Galactic Latitude/$deg$', fontsize=10)
plt.ylabel(r'Flux/ph cm^-2 s^-1', fontsize=10)
plt.xlim([-5, 5])
plt.tick_params(axis='x', labelsize=10)
plt.grid(b=True, which='major', color='0.75', linewidth=0.5)
plt.legend(prop={'size':8}, loc='lower right')
filename=raw_input('Output filename: ')
plt.savefig(filename)
