"""Functions for computing source images
from catalogs of sources with given parameters.

@todo: use xml catalogs like Fermi to be more flexible"""
import logging
import numpy as np
from gammapy.image import utils

__all__ = ['to_image']


def _add_source(image, morph_type, pars, l, b):
    """Adds a specified source at a given position."""
    from gammapy.morphology.shapes import morph_types
    # Get position and extension info
    # TODO: this won't work for delta2d, which doesn't
    # have a pars[3] entry!
    glon, glat, ext = pars[0], pars[1], pars[3]
    logging.debug('Adding source of type {0} at position '
                 '{1}, {2} with extension {3}'
                 ''.format(morph_type, glon, glat, ext))

    # Get indices of position
    #import IPython; IPython.embed()
    x_pix, y_pix = np.arange(image.header['NAXIS1']), np.arange(image.header['NAXIS2'])

    # Correct pixel numbering
    x_pix = np.floor(x_pix) - 1
    y_pix = np.floor(y_pix) - 1

    # Get data
    data = image.data
    pixsize = image.header['CDELT2']

    # Compute boxsize with one pixel margin
    boxsize = 2 * np.floor(ext / pixsize + 2)

    # Determine corner pixels
    left_bottom = [y_pix, x_pix] - boxsize / 2
    right_top = [y_pix, x_pix] + boxsize / 2

    # Determine image boundaries
    y_bound, x_bound = data.shape

    # Correct Behaviour at the boundaries
    if left_bottom[1] < 0 or left_bottom[0] < 0 or right_top[1] > x_bound or right_top[0] > y_bound:
        logging.debug('Source is not completely in image.')

        # Get out of range coordinate map
        l, b = utils.coordinates(image, x_pix, y_pix, boxsize)

        # Get image of source
        box_image = morph_types[morph_type](pars, l, b)

        # Get slices
        s_x, s_y, b_x, b_y = _get_slices(left_bottom, right_top, x_bound, y_bound)

        # For debugging
        logging.debug('left_bottom: {0}'.format(left_bottom))
        logging.debug('right_top: {0}'.format(right_top))
        logging.debug('box_image.shape: {0}'.format(box_image.shape))
        logging.debug('s_x: {0}'.format(s_x))
        logging.debug('s_y: {0}'.format(s_y))
        logging.debug('b_x: {0}'.format(b_x))
        logging.debug('b_y: {0}'.format(b_y))
        logging.debug('box_image.shape: {0}'.format(box_image[b_y, b_x].shape))
        logging.debug('data.shape: {0}'.format(data[s_y, s_x].shape))

        # Add box to image
        data[s_y, s_x] += box_image[b_y, b_x]

    else:
        logging.debug('Source is completely in image.')
        s_x = slice(left_bottom[1], right_top[1])
        s_y = slice(left_bottom[0], right_top[0])
        # Get image of source
        box_image = morph_types[morph_type](pars, l[s_y, s_x], b[s_y, s_x])
        # Add box to image
        data[s_y, s_x] += box_image


def _get_slices(left_bottom, right_top, x_bound, y_bound):
    # Set up slices in x direction
    if left_bottom[1] < 0:
        logging.debug('Left side out of range')
        s_x = slice(0, right_top[1])
        b_x = slice(-left_bottom[1], right_top[1] - left_bottom[1])

    elif right_top[1] > x_bound:
        logging.debug('Right side out of range')
        s_x = slice(left_bottom[1], x_bound)
        b_x = slice(0, x_bound - left_bottom[1])

    elif left_bottom[1] < 0 and right_top[1] > x_bound:
        s_x = slice(0, x_bound)
        b_x = slice(-left_bottom[1], x_bound - left_bottom[1])

    else:
        s_x = slice(left_bottom[1], right_top[1])
        b_x = slice(0, right_top[1] - left_bottom[1])

    # Set up slices in y direction
    if left_bottom[0] < 0:
        logging.debug('Bottom side out of range')
        s_y = slice(0, right_top[0])
        b_y = slice(-left_bottom[0], right_top[0] - left_bottom[0])

    elif right_top[0] > y_bound:
        logging.debug('Top side out of range')
        s_y = slice(left_bottom[0], y_bound)
        b_y = slice(0, y_bound - left_bottom[0])

    elif left_bottom[0] < 0 and right_top[0] > y_bound:
        s_y = slice(0, y_bound)
        b_y = slice(-left_bottom, y_bound - left_bottom[0])

    else:
        s_y = slice(left_bottom[0], right_top[0])
        b_y = slice(0, right_top[0] - left_bottom[0])

    return s_x, s_y, b_x, b_y


def _to_image_simple(catalog, image):
    """Add sources from a catalog to an image.

    catalog = atpy.Table
    image = kapteyn.maputils.FITSimage

    Note: This implementation doesn't use bounding boxes and thus is slow.
    You should use the faster _to_image_bbox()"""
    from kapteyn.maputils import FITSimage
    from morphology.shapes import morph_types
    nsources = len(catalog)
    data = image.dat
    logging.info('Adding {0} sources.'.format(nsources))

    # Get coordinate maps
    l, b = utils.coordinates(image)

    # Add sources to image one at a time
    for i in range(nsources):
        logging.debug('Adding source {0:3d} of {1:3d}.'
                      ''.format(i, nsources))
        # nans = np.isnan(data).sum()
        # if nans:
        #    logging.debug('Image contains {0} nan entries.'.format(nans))

        # Get relevant columns
        morph_type = catalog[i]['morph_type']
        xpos = catalog[i]['glon']
        xpos = np.where(xpos > 180, xpos - 360, xpos)
        ypos = catalog[i]['glat']
        ampl = catalog[i]['ampl']
        sigma = catalog[i]['sigma']
        epsilon = catalog[i]['epsilon']
        theta = catalog[i]['theta']
        r_in = catalog[i]['r_in']
        r_out = catalog[i]['r_out']

        pars = {'delta2d': (xpos, ypos, ampl),
                'shell2d': (xpos, ypos, ampl, r_in, r_out),
                'gauss2d': (xpos, ypos, ampl, sigma, epsilon, theta)}
        par = pars[morph_type]

        data += morph_types[morph_type](par, l, b)

    return FITSimage(externalheader=image.hdr, externaldata=data)


def _to_image_bbox(catalog, image):
    """@todo: make the usage of bbox an option of the function
    _to_image_simple, otherwise we duplicate a lot of code!


    Function takes catalog and an empty image.
    It returns the image of the catalog.
    catalog = atpy.Table
    image = kapteyn.FITSimage

    @todo: This function should have NO explicit references to
    morphology column names.
    Do something like:
    for p in morphology.morph_pars:
        exec('p = {0}'.format(p))

    and then
    pars = [eval('[??? for x in eval('morphology.{0}_par'.format(???
    """
    logging.info('Building image out of catalog.')
    # Select sources in FOV
    # TODO: improve this check by reading the min, max from
    # the FITS header instead of finding the min of the array.
    l, b = utils.coordinates(image)#, glon_sym=True)

    # Get catalog columns common to all morph types
    morph_type = catalog.field('morph_type')
    glon = catalog.field('glon')
    glat = catalog.field('glat')
    S_PWN = catalog.field('S_PWN')
    S_SNR = catalog.field('S_SNR')
    ext_in_SNR = catalog.field('ext_in_SNR')
    ext_out_SNR = catalog.field('ext_out_SNR')
    ext_out_PWN = catalog.field('ext_out_PWN')

    contributing = np.all([l.min() < glon + ext_out_SNR, glon - ext_out_SNR < l.max(),
                           b.min() < glat + ext_out_SNR, glat - ext_out_SNR < b.max(), S_SNR > 0], axis=0)

    # Make image, adding one source at a time
    nsources = contributing.sum()
    logging.info('Adding {0} sources.'.format(nsources))
    logging.info('Skipping {0} sources.'.format(catalog.size - nsources))
    for i in range(catalog.size):
        # Skip sources that are not contributing
        if not contributing[i]:
            # logging.debug('Skipping source {0:6d}.'.format(i))
            continue

        # Make a list of parameters appropriate for the morph_type
        pars = [glon[i], glat[i]]
        if morph_type[i] == 'sphere2d':
            pars += [S_PWN[i], ext_out_PWN[i]]
        elif morph_type[i] == 'shell2d':
            pars += [S_SNR[i], ext_out_SNR[i], ext_in_SNR[i]]

        logging.debug('Adding source {0:6d}.'.format(i))
        _add_source(image, morph_type[i], pars, l, b)


def _to_image_gaussian(catalog, image):
    """Add catalog of asymmetric Gaussian sources to an image.

    @type catalog: atpy.Table
    @type image: maputils.FITSimage
    @return: maputils.FITSimage"""
    from kapteyn.maputils import FITSimage
    logging.info('Adding %d sources.' % len(catalog.data))

    #
    # Read catalog
    #
    l_center = catalog.data['glon']
    b_center = catalog.data['glat']
    a = catalog.data['a']
    b = catalog.data['b']
    theta = catalog.data['theta']
    flux = catalog.data['flux']

    # Compute matrix entries that describe the ellipse, as defined in
    # the section "9.1.6. Ellipse parameters" in the SExtractor Manual.
    # Note that clb already contains the factor of 2!
    #
    # Another reference giving formulas is
    # http://en.wikipedia.org/wiki/Gaussian_function
    cll = np.cos(theta) ** 2 / a ** 2 + np.sin(theta) ** 2 / b ** 2
    cbb = np.sin(theta) ** 2 / a ** 2 + np.cos(theta) ** 2 / b ** 2
    clb = 2 * np.cos(theta) * np.sin(theta) * (a ** -2 - b ** -2)

    #
    # Add sources to image
    #
    data = image.dat
    header = image.hdr
    for i in range(len(catalog)):  # Loop sources
        l, b = utils.coordinates(image)
        # convert l to the range -180 to +180
        l = np.where(l > 180, l - 360, l)
        exponent = cll[i] * (l - l_center[i]) ** 2 + \
            cbb[i] * (b - b_center[i]) ** 2 + \
            clb[i] * (l - l_center[i]) * (b - b_center[i])
        data += flux[i] * np.exp(-exponent)

    # Return image with sources
    image_out = FITSimage(externalheader=header,
                          externaldata=data)
    return image_out

# Select which to_image method is the default here:
to_image = _to_image_simple
# to_image = _to_image_bbox
# to_image = _to_image_gaussian
