import os
import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval
import matplotlib.pyplot as plt
from matplotlib import rcParams


def plot_images(image1, image2, group=0, title1=None, title2=None):
    """Display two images side-by-side for comparison. 

    As a quick-look function, this is not robust or versatile. For
    example, it uses ZScaleInterval for all image normalization. It
    shows the first integration for all images, appropriate for CEERS
    exposures that only have 1 integration.

    Args:
        image1 (str): filename of first image for comparison
        image2 (str): filename of second image for comparison
        group (Optional [int]): index of group number to plot for the
            4D uncal images
        title1 (Optional [str]): title for image1 plot
        title2 (Optional [str]): title for image2 plot
    """
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 10
    rcParams['axes.titlesize'] = 25


    im1 = fits.getdata(image1, 'SCI')
    im2 = fits.getdata(image2, 'SCI')

    # If images are 4D, pick slices to plot
    if len(im1.shape) == 4:
        # take first integration - there's only 1 for CEERS images
        im1 = im1[0,group,:,:]
    if len(im2.shape) == 4:
        im2 = im2[0,group,:,:]

    # normalize with zscale
    norm1 = ImageNormalize(im1, interval=ZScaleInterval())
    norm2 = ImageNormalize(im2, interval=ZScaleInterval())

    fig,(ax1,ax2) = plt.subplots(1, 2, tight_layout=True)
    ax1.imshow(im1, origin='lower', interpolation='none', cmap='Greys',
               norm=norm1)
    ax2.imshow(im2, origin='lower', interpolation='none', cmap='Greys',
               norm=norm2)

    for ax in [ax1,ax2]:
        ax.set_xlabel('pixels', fontsize=10)
    ax1.set_ylabel('pixels', fontsize=10)


    if title1:
        ax1.set_title(title1)
    if title2:
        ax2.set_title(title2)


