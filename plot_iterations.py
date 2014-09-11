"""Script to plot background and significance images for each iteration of the source/background separation code"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

n_iterations = 5

#plt.clf()
plt.figure(figsize=(20, 5))

for iteration in range(n_iterations):
    filename = 'test{0:02d}_mask_true_region.fits'.format(iteration)
    mask = fits.getdata(filename)#[100:300,:]

    plt.subplot(n_iterations, 2, 2 * iteration + 1)
    filename = 'test{0:02d}_background_true_region.fits'.format(iteration)
    data = fits.getdata(filename)#[100:300,:]
    plt.imshow(data)
    plt.contour(mask, levels=[0], linewidths=2, colors='white')
    plt.axis('off')
    plt.tight_layout()
    #plt.title(filename)
    
    plt.subplot(n_iterations, 2, 2 * iteration + 2)
    filename = 'test{0:02d}_significance_true_region.fits'.format(iteration)
    data = fits.getdata(filename)#[100:300,:]
    plt.imshow(data, vmin=-3, vmax=5)
    plt.contour(mask, levels=[0], linewidths=2, colors='white')
    plt.axis('off')
    plt.tight_layout()
    #plt.title(filename)
    #plt.colorbar()

plt.tight_layout()
plt.savefig('iterations.pdf')