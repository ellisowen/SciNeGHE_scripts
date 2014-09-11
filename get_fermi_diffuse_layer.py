""" Gets correct log cneter energy bin image for the fermi diffuse model"""

from scipy.stats import gmean
import numpy as np
from astropy.io import fits

fd_cube = fits.open('gll_iem_v05.fits')
energies = fd_cube[1].data

energy = gmean([10, 500]) * 1000 #GeV to MeV

energy_array = fd_cube[1].data['Energy']
index = np.searchsorted(energy_array, energy)

actual_energy = energy_array[index]
print actual_energy
from kapteyn.maputils import FITSimage
from gammapy.image.utils import cube_to_image
energy_image = cube_to_image(fd_cube[0], index)
energy_image.writeto('fermi_diffuse_slice_10_500.fits', clobber=True)
#reproject to 0.1 deg per pixel resolution
original = FITSimage('fermi_diffuse_slice_10_500.fits', hdunr=1)
reference_data = fits.open('10_500_counts.fits')[0].data
reference_header = fits.open('10_500_counts.fits')[0].header
reference = fits.ImageHDU(data = reference_data, header = reference_header)
reference.writeto('10_500_counts2.fits', clobber=True)
reference = FITSimage('10_500_counts2.fits', hdunr=1)

from utils import image_reprojection

new = image_reprojection(original, reference)
new.writetofits('fermi_diffuse_10_500_reprojected.fits', clobber=True)
out_data = np.nan_to_num(fits.open('fermi_diffuse_10_500_reprojected.fits')[0].data)
out_header = fits.open('fermi_diffuse_10_500_reprojected.fits')[0].header
out = fits.ImageHDU(data = out_data, header = out_header)
out.writeto('fermi_diffuse_10_500_reprojected.fits', clobber=True)


