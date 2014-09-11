"""Estimate a diffuse emission model from Fermi LAT data.
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.units import Quantity
from gammapy.datasets import FermiGalacticCenter
from gammapy.background import IterativeKernelBackgroundEstimator, GammaImages
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image import make_empty_image, catalog_image, binary_disk
from gammapy.image.utils import cube_to_image, solid_angle
from gammapy.data import SpectralCube
from gammapy.image.utils import WCS
from gammapy.spectrum.flux_point import _energy_lafferty_power_law

# *** PREPARATION ***

# Parameters

CORRELATION_RADIUS = 3 # pix
SIGNIFICANCE_THRESHOLD = 4
MASK_DILATION_RADIUS = 3

psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)

# *** LOADING INPUT ***

# Counts must be provided as a counts ImageHDU
flux_file = raw_input('Counts Map: ')
exposure_file = raw_input('Exposure Map: ')
spec_ind = input('Spectral Index (for reprojection): ')
flux_hdu = fits.open(flux_file)[0]
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

counts_data = flux_spec_cube.data # Input is already a counts cube in this case
counts = fits.ImageHDU(data=counts_data[0], header=flux_hdu.header)
# Start with flat background estimate
# Background must be provided as an ImageHDU
background_data=np.ones_like(counts_data, dtype=float)
background = fits.ImageHDU(data=background_data[0], header=flux_hdu.header)
images = GammaImages(counts=counts, background=background)

source_kernel = binary_disk(CORRELATION_RADIUS).astype(float)
source_kernel /= source_kernel.sum()

background_kernel = np.ones((5, 100))
background_kernel /= background_kernel.sum()

# *** ITERATOR ***

ibe = IterativeKernelBackgroundEstimator(images=images,
                                         source_kernel=source_kernel,
                                         background_kernel=background_kernel,
                                         significance_threshold=SIGNIFICANCE_THRESHOLD,
                                         mask_dilation_radius=MASK_DILATION_RADIUS
                                         )

mask, new_background = ibe.run()

flux_background_data = new_background.data / exposure_spec_cube.data[0]

flux_background = fits.ImageHDU(data=flux_background_data.value, header=flux_hdu.header)

filebase = raw_input('Output file base: ')

flux_background.writeto('{0}_background.fits'.format(filebase), clobber=True)
mask.writeto('{0}_mask.fits'.format(filebase), clobber=True)

