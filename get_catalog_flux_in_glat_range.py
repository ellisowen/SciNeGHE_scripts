"""Script to compute the flux in sources from a catalog within a given latitude and longitude range"""

from astropy.io import fits
from gammapy.datasets import fetch_fermi_catalog


#filename = 'simulated_galaxy_3.fits'
glat_min = -5
glat_max = +5
glon_min = -100
glon_max = +100


#data = fits.open(filename)[1].data
data = fetch_fermi_catalog('1FHL', 'LAT_Point_Source_Catalog')
glats = data['GLAT']
glons = data['GLON']
fluxes = data['Flux']

flux_threshold = 5e-10

mask_lat = (glat_min < glats) & (glats < glat_max)
mask_lon = (glon_min < glons) & (glons < glon_max) 
mask_total = mask_lat & mask_lon
mask_flux = (flux_threshold<fluxes) & mask_total





total_flux = data['Flux'][mask_total].sum()

print total_flux