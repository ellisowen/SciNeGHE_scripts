from __future__ import print_function, division
import numpy as np

__all__ = ['mc']


def mc(model, xlim, yerr, npoints=5, verbose=0):
    """Generate points with equal-log spacing in x given by the model.

    model = [function, parameters, constants]
    xlim = [xmin, xmax]
    yerr = [ydn_err, yup_err] = fractional error on model y
    Returns:
    data = x, y, ydn, yup
    """
    # Generate equal-log spaced x points
    logx = np.linspace(np.log10(xlim[0]), np.log10(xlim[1]), npoints)
    x = 10 ** logx

    # Unpack model components
    f, p, c = model

    # Compute true y and asymmetric errors
    y = f(p, c)
    ydn = y * yerr[0]
    yup = y * yerr[1]

    #
    # Compute observed y by drawing a random value
    #

    # First decide if an up or down fluctuation occurs
    fluctuate_up = np.random.randint(0, 2, npoints)  # 1 = yes, 2 = no

    # Then draw a random y value
    yobs_dn = np.fabs(np.random.normal(0, ydn, size=npoints))
    yobs_up = np.fabs(np.random.normal(0, yup, size=npoints))
    yobs = y + np.where(fluctuate_up == 1, yobs_up, -yobs_dn)

    if verbose > 0:
        for i in range(npoints):
            fmt = '%2d' + ' %10g' * 7
            vals = i, x[i], y[i], ydn[i], yup[i], yobs_dn[i], yobs_up[i], yobs[i]
            print(fmt.format(vals))
    data = x, yobs, ydn, yup
    return data