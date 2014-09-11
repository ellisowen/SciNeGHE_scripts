"""Estimate a diffuse emission model from Fermi LAT data.
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.utils.data import download_file
from gammapy.datasets import FermiGalacticCenter
from gammapy.background import IterativeKernelBackgroundEstimator, GammaImages
from gammapy.irf import EnergyDependentTablePSF
from gammapy.image import make_empty_image, catalog_image, binary_disk
from gammapy.image.utils import cube_to_image, solid_angle

counts = fits.open('fermi_counts_galactic.fits')[0]

# Start with flat background estimate
# Background must be provided as an ImageHDU
background_data=np.ones_like(counts.data, dtype=float)
background = fits.ImageHDU(data=background_data, header=counts.header)
images = GammaImages(counts=counts, background=background)

source_kernel = binary_disk(3).astype(float)
source_kernel /= source_kernel.sum()

background_kernel = np.ones((5, 50))
background_kernel /= background_kernel.sum()

# *** ITERATOR ***

ibe = IterativeKernelBackgroundEstimator(images=images,
                                         source_kernel=source_kernel,
                                         background_kernel=background_kernel,
                                         significance_threshold=4,
                                         mask_dilation_radius=3,
                                         save_intermediate_results=True
                                         )

ibe.run(filebase='TestOutput')
#n_iterations = 3
"""
# *** RUN & PLOT ***
plt.figure(figsize=(4,6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.subplot(n_iterations+3, 1, 1)
background_hdu = counts
data = counts.data[:, 700:1400]
plt.imshow(data, vmin=-3, vmax=5)
plt.axis('off')
plt.title(r'Input Counts Image', fontsize='small')

ibe.run_iteration()

plt.subplot(n_iterations+3, 1, 2)
background_hdu = ibe.background_image_hdu
data = background_hdu.data[:, 700:1400]
plt.imshow(data)#, vmin=0, vmax=0.2)
plt.axis('off')
plt.title(r'Initial Background Estimation', fontsize='small')

for iteration in range(1, n_iterations):

    mask_hdu = ibe.mask_image_hdu
    mask = mask_hdu.data[:, 700:1400]
    mask_hdu.writeto('mask_{0}.fits'.format(iteration))
    
   # import IPython; IPython.embed()

    plt.subplot(n_iterations+3, 1, 2*iteration+1)
    significance_hdu = ibe.significance_image_hdu
    data = significance_hdu.data[:, 700:1400]
    plt.imshow(data, vmin=-3, vmax=5)
    plt.contour(mask, levels=[0], linewidths=1, colors='black')
    plt.axis('off')
    plt.title(r'Significance Image, Iteration {0}'.format(iteration), fontsize='small')

    ibe.run_iteration()

    plt.subplot(n_iterations+3, 1, 2*iteration+2)
    background_hdu = ibe.background_image_hdu
    data = background_hdu.data[:, 700:1400]
    plt.imshow(data)#, vmin=0, vmax=0.2)
    plt.contour(mask, levels=[0], linewidths=1, colors='white')
    plt.axis('off')
    plt.title(r'Background Estimation, Iteration {0}'.format(iteration), fontsize='small')

plt.tight_layout()
plt.savefig('iterations.pdf')
"""
