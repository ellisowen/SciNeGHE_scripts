"""Produces latitude profile for effectiveness of reconstruction, source emssion and background emission in latitude bins"""

from utils import compute_longitude_profile
from class_source_diffuse_estimation import RunAlgorithm, source_image, diffuse_image, data_image, correlation_radius, significance_threshold, mask_dilation_radius, LonHigh, LonLow
from grid_binning import GridSeparation
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

__all__ = ['Run', 'BackgroundProfile', 'SignificanceProfile']

EXPOSURE = 1.5768e+11

from gammapy.image.utils import solid_angle
reference = fits.open(data_image)
PIXEL_SA = solid_angle(reference).mean()

def Run(correlation_radius = correlation_radius, significance_threshold=significance_threshold, mask_dilation_radius=mask_dilation_radius, number_of_iterations=50, data_image = data_image):
    
    operate = RunAlgorithm(analysis='Data', correlation_radius=correlation_radius, significance_threshold=significance_threshold, mask_dilation_radius=mask_dilation_radius,
                           number_of_iterations=number_of_iterations, data_image=data_image)
    
    background = fits.open(operate.background_filename)
    mask = fits.open(operate.mask_filename)
    significance = fits.open(operate.significance_filename)
    
    return background, mask, significance


def BackgroundProfile(binsz):
    background, mask, signifiance = Run(correlation_radius = correlation_radius, significance_threshold=significance_threshold, mask_dilation_radius=mask_dilation_radius, number_of_iterations=50, data_image = data_image)
    background = background[1]
    glons = np.arange(LonLow, LonHigh+binsz, binsz)    
    glon_bounds = Table()
    glon_bounds['CHANNEL'] = np.arange(len(glons) - 1)
    glon_bounds['GLON_MIN'] = np.float64(glons[:-1])
    glon_bounds['GLON_MAX'] = np.float64(glons[1:])
    data = compute_longitude_profile(glon_bounds = glon_bounds, binsz=binsz, image=background, datatype=2, emission=4, tev=0)
    return data

def GeneralProfile(binsz, filename):
    """This should be improved to avoid hard-coded filenames"""
    data_image = fits.open(filename)
    data_image= data_image[1]
    glons = np.arange(LonLow, LonHigh+binsz, binsz)    
    glon_bounds = Table()
    glon_bounds['CHANNEL'] = np.arange(len(glons) - 1)
    glon_bounds['GLON_MIN'] = np.float64(glons[:-1])
    glon_bounds['GLON_MAX'] = np.float64(glons[1:])
    data = compute_longitude_profile(glon_bounds = glon_bounds, binsz=binsz, image=data_image, datatype=2, emission=4, tev=0)
    return data

def PlotProfile(label1, label2, data1, data2, filename, smooth, binsz):
    plt.clf()
    plt.figure(figsize=(10,4))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    
    
    
    counts1 = np.append(data1['FLUX'].data/(EXPOSURE*PIXEL_SA), 0)
    glons1 = np.append(data1['GLON_MIN'][0], data1['GLON_MAX'].data)

    profile1 = np.histogram(glons1, bins=np.sort(glons1), weights=counts1)
    xstep=binsz
    y_smooth_1 = gaussian_filter(profile1[0], smooth / xstep)
    print y_smooth_1.mean()
    
    counts2 = np.append(data2['FLUX'].data/(EXPOSURE*PIXEL_SA), 0)
    glons2 = np.append(data2['GLON_MIN'][0], data2['GLON_MAX'].data)

    profile2 = np.histogram(glons2, bins=np.sort(glons2), weights=counts2)
    xstep=binsz
    y_smooth_2 = gaussian_filter(profile2[0], smooth / xstep)
    print y_smooth_2.mean()
    
    x1 = 0.5 * (glons1[1:] + glons1[:-1])
    if abs(LonHigh) != abs(LonLow): #TODO: Fix this hack!
        plt.plot((-1*x1)-15, y_smooth_1, label='{0}'.format(label1))
        plt.plot((-1*x1)-15, y_smooth_2, label='{0}'.format(label2))
    else:
        plt.plot((-1*x1), y_smooth_1, label='{0}'.format(label1))
        plt.plot((-1*x1), y_smooth_2, label='{0}'.format(label2))
    plt.xlabel(r'Galactic Longitude/$deg$', fontsize=10)
    plt.ylabel(r'Surface Brightness/ph cm$^{-2}$ s$^{-1} sr^{-1}$', fontsize=10)
    plt.xlim([LonHigh, LonLow])
    plt.tick_params(axis='x', labelsize=10)
    plt.grid(b=True, which='major', color='0.75', linewidth=0.5)
    plt.legend(prop={'size':8})
    plt.savefig(filename)
    
#def PlotProfile2(label, data, filename, smooth, binsz):
#    plt.clf()
#    plt.figure(figsize=(5,4))
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    counts = np.append(data['FLUX'].data/EXPOSURE, 0)
#    glons = np.append(data['GLON_MIN'][0], data['GLON_MAX'].data)
#    profile = np.histogram(glons, bins=np.sort(glons), weights=counts)
#    xstep=binsz
#    y_smooth_1 = gaussian_filter(profile[0], smooth / xstep)
#    x1 = 0.5 * (glons[1:] + glons[:-1])
#    plt.plot(x1, y_smooth_1, label='{0}'.format(label))
#    plt.hlines(0, -1, 1)
#    plt.xlabel(r'Galactic Longitude/$deg$', fontsize=10)
#    plt.ylabel(r'Relative Residual Significance', fontsize=10)
#    plt.xlim([-1, 1])
#    plt.tick_params(axis='x', labelsize=10)
#    plt.tick_params(axis='y', labelsize=10, colors='w')
#    plt.grid(b=True, which='major', color='0.75', linewidth=0.5)
#    plt.savefig(filename)
    
    
#def SignificanceProfile(binsz):
#    background, mask, significance = Run(correlation_radius = correlation_radius, significance_threshold=significance_threshold, mask_dilation_radius=mask_dilation_radius, number_of_iterations=50, data_image = data_image)
#    significance = significance[1]
#    glons = np.arange(-1.1, 1.1+binsz, binsz)    
#    glon_bounds = Table()
#    glon_bounds['CHANNEL'] = np.arange(len(glons) - 1)
#    glon_bounds['GLON_MIN'] = np.float64(glons[:-1])
#    glon_bounds['GLON_MAX'] = np.float64(glons[1:])
#    data = compute_latitude_profile(glon_bounds = glon_bounds, binsz=binsz, image=significance, datatype=2, emission=4, tev=0)
#    return data

if __name__ == '__main__':
    
    binsz=0.1
    background_data = BackgroundProfile(binsz)
    initial_image_data = GeneralProfile(binsz, 'data_image_true_size.fits')
    initial_background_data = GeneralProfile(binsz, 'background_input_true_size.fits') 
    PlotProfile('Estimated Background', 'Initial Total Image', background_data, initial_image_data, 'BackgroundLONProfile_small_data.pdf', 0.2, binsz)
    
