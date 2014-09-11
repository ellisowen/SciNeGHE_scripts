""" Script to determine what fraction of source and diffuse flux is found in the large and small region"""

from astropy.io import fits

filename = raw_input('Image filename: ')
source = fits.open(filename)[1]

from gammapy.image.utils import coordinates, solid_angle
lats, lons = coordinates(source)
solid_angle = solid_angle(source)
region_glat_high = 5
region_glat_low = -5
region_glon_high = 100
region_glon_low = -100

mask1 = (region_glat_low < lats) & (lats < region_glat_high)
mask2 = (region_glon_low < lons) & (lons < region_glon_high)
mask = mask1 & mask2

a = source.data #/ solid_angle.value
source_flux_frac = a[mask].sum()/a.sum()
print "Total Flux in Image"
print a.sum()

print "Galactic Flux"
print a[mask].sum()

print "Source Flux Fraction"
print source_flux_frac

