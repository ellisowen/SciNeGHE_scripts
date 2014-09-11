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
SIGNIFICANCE_THRESHOLD = 5
MASK_DILATION_RADIUS = 0.3

psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)

# *** LOADING INPUT ***

# Counts must be provided as a counts ImageHDU
fermi_diffuse = 'gll_iem_v05_rev1.fit'
comparison_diffuse = raw_input('Diffuse Model: ')

spec_ind = input('Spectral Index (for reprojection): ')
flux_hdu = fits.open(comparison_diffuse)[1]
flux_wcs = WCS(flux_hdu.header)
energy_flux = Quantity([_energy_lafferty_power_law(10000, 500000, spec_ind)], 'MeV')
flux_data = np.zeros((1, 103, 2001))
flux_data[0] = Quantity(flux_hdu.data, '')
flux_spec_cube = SpectralCube(data=Quantity(flux_data, '1 / (cm2 MeV s sr)'), wcs=flux_wcs, energy=energy_flux)

diffuse_spec_cube = SpectralCube.read(fermi_diffuse)
diffuse_spec_cube.data = Quantity(diffuse_spec_cube.data, '1 / (cm2 MeV s sr)')
fermi_diffuse_spec_cube = diffuse_spec_cube.reproject_to(flux_spec_cube)

energy_band = Quantity((10, 500), 'GeV')

a = fermi_diffuse_spec_cube.integral_flux_image(energy_band)
angles = fermi_diffuse_spec_cube.solid_angle_image
b = angles.value
c = a.data * b

diffuse = fits.ImageHDU(data=c, header=flux_hdu.header)
diffuse.writeto('diffuse_correct.fits')

ratio = flux_spec_cube.data / c

ratio_hdu = fits.ImageHDU(data=ratio.value, header=flux_hdu.header)
ratio_im_hdu = cube_to_image(ratio_hdu)
ratio_im_hdu.writeto('ratio.fits', clobber=True)
