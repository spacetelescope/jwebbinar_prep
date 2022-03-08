import os
import shutil
from glob import glob
import argparse
import logging
from datetime import datetime
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
import yaml
# jwst-related imports
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.flatfield.flat_field import do_correction
from stdatamodels import util
import crds
# logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


"""
Measure striping by collapsing image over rows and columns, using the
sigma-clipped median value to avoid source flux. 

The measurement/subtraction needs to be done along one axis at a time,
since the measurement along x will depend on what has been subtracted
from y.
"""


def gaussian(x, a, mu, sig):
    return a * np.exp(-(x-mu)**2/(2*sig**2))


def fit_sky(data):
    """Fit distribution of sky fluxes with a Gaussian"""
    bins = np.arange(-0.5, 0.5, 0.001)
    h,b = np.histogram(data, bins=bins)
    bc = 0.5 * (b[1:] + b[:-1])
    binsize = b[1] - b[0]

    p0 = [10, bc[np.argmax(h)], 0.01]
    popt,pcov = curve_fit(gaussian, bc, h, p0=p0)

    return popt[1]


def collapse_image(im, mask, dimension='y', sig=2.):
    """collapse an image along one dimension to check for striping.

    By default, collapse columns to show horizontal striping (collapsing
    along columns). Switch to vertical striping (collapsing along rows)
    with dimension='x' 

    Striping is measured as a sigma-clipped median of all unmasked pixels 
    in the row or column.

    Args:
        im (float array): image data array
        mask (bool array): image mask array, True where pixels should be 
            masked from the fit (where DQ>0, source flux has been masked, etc.)
        dimension (Optional [str]): specifies which dimension along which 
            to collapse the image. If 'y', collapses along columns to 
            measure horizontal striping. If 'x', collapses along rows to 
            measure vertical striping. Default is 'y'
        sig (Optional [float]): sigma to use in sigma clipping
    """
    # axis=1 results in array along y
    # axis=0 results in array along x

    if dimension == 'y':
#        collapsed = np.median(im, axis=1)
        res = sigma_clipped_stats(im, mask=mask, sigma=sig, cenfunc='median',
                                        stdfunc='std', axis=1)

    elif dimension == 'x':
#        collapsed = np.median(im, axis=0)
        res = sigma_clipped_stats(im, mask=mask, sigma=sig, cenfunc='median',
                                        stdfunc='std', axis=0)

    return res[1]
    

def measure_striping(image, apply_flat=True, mask_sources=True, seedim_directory='./', threshold=0.01):
    """Removes striping in rate.fits files before flat fielding.

    Measures and subtracts the horizontal & vertical striping present in 
    countrate images. The striping is most likely due to 1/f noise, and 
    the RefPixStep with odd_even_columns=True and use_side_ref_pixels=True
    does not fully remove the pattern, no matter what value is chosen for 
    side_smoothing_length. There is also residual vertical striping in NIRCam 
    images simulated with Mirage.

    Note: 
        The original rate image file is copied to *_rate_orig.fits, and 
        the rate image with the striping patterns removed is saved to 
        *_rate.fits, overwriting the input filename

    Args:
        image (str): image filename, including full relative path
        apply_flat (Optional [bool]): if True, identifies and applies the 
            corresponding flat field before measuring striping pattern. 
            Applying the flat first allows for a cleaner measure of the 
            striping, especially for the long wavelength detectors. 
            Default is True.
        mask_sources (Optional [bool]): If True, masks out sources in image
            before measuring the striping pattern so that source flux is 
            not included in the calculation of the sigma-clipped median.
            Sources are identified using the Mirage seed images.
            Default is True.
        seedim_directory (Optional [bool]): Directory containing 
            Mirage seed images, used if mask_sources is True. 
            Default is working directory.
        threshold (Optional [float]): threshold (in ADU/s) to use in the 
            seed images when identifying pixels to mask. This will depend on 
            the seed image and brightness of input sources. Default is 0.01
    """
    try:
        crds_context = os.environ['CRDS_CONTEXT']
    except KeyError:
        crds_context = crds.get_default_context()

    model = ImageModel(image)
    log.info('Measuring image striping')
    log.info('Working on %s'%image)

    # check that striping hasn't already been removed
    for entry in model.history:
        for k,v in entry.items():
            if 'Removed horizontal,vertical striping; remstriping.py' in v:
                print('%s already cleaned. Skipping!'%os.path.basename(image))
                return

    # apply the flat to get a cleaner meausurement of the striping
    if apply_flat:
        log.info('Applying flat for cleaner measurement of striping patterns')
        # pull flat from CRDS using the current context
        crds_dict = {'INSTRUME':'NIRCAM', 
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'PUPIL':model.meta.instrument.pupil, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
        flats = crds.getreferences(crds_dict, reftypes=['flat'], 
                                   context=crds_context)
        # if the CRDS loopup fails, should return a CrdsLookupError, but 
        # just in case:
        try:
            flatfile = flats['flat']
        except KeyError:
            log.error('Flat was not found in CRDS with the parameters: {}'.format(crds_dict))
            exit()

        log.info('Using flat: %s'%(os.path.basename(flatfile)))
        with FlatModel(flatfile) as flat:
            # use the JWST Calibration Pipeline flat fielding Step 
            model,applied_flat = do_correction(model, flat)
            
    # construct mask for median calculation
    mask = np.zeros(model.data.shape, dtype=bool)
    mask[model.dq > 0] = True
    
    # mask out sources
    if mask_sources:
        log.info('Using the input seed image to mask out source flux')
        base = os.path.basename(image).split('_rate.fits')[0]
        seedimage = glob(os.path.join(seedim_directory, 
                         '%s_*_galaxy_seed_image.fits'%base))
        log.info('Using seedim: %s'%(os.path.basename(seedimage[0])))
        seedim = fits.getdata(seedimage[0])
        log.info('Masking flux above threshold %f'%threshold)
        wobj = np.where(seedim > threshold)
        mask[wobj] = True

    # measure the pedestal in the unmasked parts of the image
    log.info('Measuring the pedestal in the image')
    pedestal_data = model.data[~mask]
    pedestal_data = pedestal_data.flatten()
    median_image = np.median(pedestal_data)
    log.info('Image median (unmasked and DQ==0): %f'%(median_image))
    try:
        pedestal = fit_sky(pedestal_data)
    except RuntimeError as e:
        log.error("Can't fit sky, using median value instead")
        pedestal = median_image
    else:
        log.info('Fit pedestal: %f'%pedestal)

    # subtract off pedestal to make it easier to fit 
    model.data -= pedestal

    # fit horizontal striping, collapsing along columns
    horizontal_striping = collapse_image(model.data, mask, dimension='y')
    # remove horizontal striping, requires taking transpose of image
    temp_image = model.data.T - horizontal_striping
    # transpose back
    temp_image2 = temp_image.T

    # fit vertical striping, collapsing along rows
    vertical_striping = collapse_image(temp_image2, mask, dimension='x')

    model.close()
    
    # copy image 
    log.info('Copying input to %s'%image.replace('rate.fits', 'rate_orig.fits'))
    shutil.copy2(image, image.replace('rate.fits', 'rate_orig.fits'))

    # remove striping from science image
    with ImageModel(image) as immodel:
        sci = immodel.data
        temp_sci = sci.T - horizontal_striping
        # transpose back
        temp_sci2 = temp_sci.T
        outsci = temp_sci2 - vertical_striping
        # replace NaNs with zeros and update DQ array
        # the image has NaNs where an entire row/column has been masked out
        # so no median could be calculated.
        # All of the NaNs on LW detectors and most of them on SW detectors
        # are the reference pixels around the image edges. But there is one
        # additional row on some SW detectors 
#        refpixflag = dqflags.pixel['REFERENCE_PIXEL']
#        wref = np.bitwise_and(immodel.dq, refpixflag)
#        outsci[np.where(wref)] = 0
        wnan = np.isnan(outsci)
        bpflag = dqflags.pixel['DO_NOT_USE']
        outsci[wnan] = 0
        immodel.dq[wnan] = np.bitwise_or(immodel.dq[wnan], bpflag)

        # write output
        immodel.data = outsci
        # add history entry
        time = datetime.now()
        stepdescription = 'Removed horizontal,vertical striping; remstriping.py %s'%time.strftime('%Y-%m-%d %H:%M:%S')
        # writing to file doesn't save the time stamp or software dictionary 
        # with the History object, but left here for completeness
        software_dict = {'name':'remstriping.py',
                         'author':'Micaela Bagley',
                         'version':'1.0',
                         'homepage':'ceers.github.io'}
        substr = util.create_history_entry(stepdescription, 
                                              software=software_dict)
        immodel.history.append(substr)
        log.info('Saving cleaned image to %s'%image)
        immodel.save(image)


def main():
    parser = argparse.ArgumentParser(description=
        'Measure and remove horizontal/vertical striping from rate images')
    parser.add_argument('--output_dir', type=str, default='calibrated',
        help='Output directory for cleaned images. Default is calibrated.')
    parser.add_argument('--runone', type=str,
        help='Filename of single file to clean. If set, overrides the runall argument')
    parser.add_argument('--runall', action='store_true',
        help='Set to run all *rate.fits images in the output_dir directory.')
    parser.add_argument('--apply_flat', action='store_true',
        help='Set to apply the flat field before measuring striping pattern.')
    parser.add_argument('--mask_sources', action='store_true',
        help='Set to mask out sources in image before measuring striping pattern.')
    parser.add_argument('--seedim_directory', type=str, default='./',
        help='Directory containing Mirage seed images, used if mask_sources is True. Default is working directory.')
    parser.add_argument('--threshold', type=float, default=0.01,
        help='Threshold (in ADU/s) to use in the seed images when identifying pixels to mask. Default is 0.01')

    args = parser.parse_args()

    if args.runone:
        image = os.path.join(args.output_dir, args.runone)
        measure_striping(image, apply_flat=args.apply_flat, 
                         mask_sources=args.mask_sources,
                         seedim_directory=args.seedim_directory, 
                         threshold=args.threshold)

    elif args.runall:
        rates = glob(os.path.join(args.output_dir, '*rate.fits'))
        rates.sort()
        for rate in rates:
            measure_striping(rate, apply_flat=args.apply_flat,
                             mask_sources=args.mask_sources,
                             seedim_directory=args.seedim_directory,
                             threshold=args.threshold)

    else:
        print('Specify either --runone with a single input filename or --runall to process all *_rate.fits files')


if __name__ == '__main__':
    main()

    
