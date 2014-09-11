#!/usr/bin/env python
"""
Given a catalog of sources, simulate a flux image. 
"""

from astropy.io import fits
from gammapy.image.utils import make_empty_image
from simulate import _to_image_bbox as to_image

catalog = fits.open('mc_catalog.fits')[1].data
image = make_empty_image(nxpix=600, nypix=600, binsz=0.02, xref=0, yref=0, dtype='float64')
to_image(catalog, image)
image.writetofits('test_image.fits', clobber=True)