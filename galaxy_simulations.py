"""Simulates a galaxy of point sources and produces an image.
"""
from astropy.units import Quantity
from astropy.io import fits
from astropy.table import Table
from aplpy import FITSFigure
from gammapy.astro import population
from gammapy.datasets import FermiGalacticCenter
from gammapy.image import make_empty_image, catalog_image
from gammapy.irf import EnergyDependentTablePSF
from gammapy.utils.random import sample_powerlaw

# Create image of defined size
reference = make_empty_image(nxpix=3600, nypix=1800, binsz=0.1)

psf_file = FermiGalacticCenter.filenames()['psf']
psf = EnergyDependentTablePSF.read(psf_file)
filename = raw_input('Input simulated table filename: ')
table = fits.open(filename)[1].data
table = Table(table)
table.meta['Energy Bins'] = Quantity([10, 500], 'GeV')
# Create image

image = catalog_image(reference, psf, catalog='simulation', source_type='point',
                  total_flux=True, sim_table=table)
hdu = image.to_fits()[0]
hdu.writeto('image_{0}'.format(filename), clobber=True)

# Plot
fig = FITSFigure(image.to_fits()[0])
fig.show_grayscale(stretch='linear', interpolation='none')
fig.add_colorbar()
fig.save('testimage.pdf')
