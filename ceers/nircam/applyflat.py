import os
import shutil
from glob import glob
import argparse
import logging
from datetime import datetime
import numpy as np
# jwst-related imports
from jwst.datamodels import ImageModel, FlatModel
from jwst.flatfield.flat_field import do_correction
from stdatamodels import util
# Pipeline logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def fix_header(image):
    """Add history entry to header reflecting apply_custom_flat step

    Adding the history entry in apply_custom_flat() did not work, the 
    history in the saved file was not updated. Reappending the entry 
    after opening the file works
    """
    im = ImageModel(image)
    # add history entry
    time = datetime.now()
    stepdescription = 'Applied custom flat for NRCA5 detector; applyflat.py %s'%time.strftime('%Y-%m-%d %H:%M:%S')
    # writing to file doesn't save the time stamp or software dictionary 
    # with the History object, but left here for completeness
    software_dict = {'name':'applyflat.py',
                     'author':'Micaela Bagley',
                     'version':'1.0',
                     'homepage':'ceers.github.io'}
    substr = util.create_history_entry(stepdescription,
                                          software=software_dict)
    new_history = []
    for item in im.history:
        new_history.append(item)
    new_history.append(substr)
    im.history = new_history
    im.save(image)


def apply_custom_flat(image, suffix='_unflat.fits', data_dir='.'):
    """Apply a custom flat to NRCA5 cal images to remove an extra feature

    The feature is likely the result of reference file flipping about x
    at some point during the Mirage image creation, and so it is not 
    properly removed by the pipeline. 

    The custom flats have been specfically created for the CEERS 5 
    pointing and observation specification. They are the result of 
    combining 30 empty input files in each filter.

    apply_custom_flat uses jwst.flatfield.flat_field.do_correction() 
    to correct the cal files with the user-specified flat field. This 
    handles checking for DQ flags and propagating the error arrays

    Note: 
        The original cal file is copied to *_[suffix].fits, and cal
        file with the feature removed is saved to [image], 
        overwriting the input filename

    Args:
        image (str): NRCA5 *cal.fits file, including full relative path
        suffix (Optional [str]): suffix to add to original input files 
    """
    im = ImageModel(image)
    
    # check that flat hasn't already been applied
    for entry in im.history:
        for k,v in entry.items():
            if 'Applied custom flat for NRCA5 detector; applyflat.py' in v:
                print('%s flat already applied. Skipping!'%os.path.basename(image))
                return


    # get image filter
    filt = im.meta.instrument.filter
    log.info('Applying custom flat for NRCA5 %s'%filt)
    customflat = os.path.join(data_dir, 'customflats/ceers5_customflat_nrca5_%s.fits'%filt.lower())
    log.info('Custom flat: %s'%customflat)
    log.info('Applying to %s'%image)

    # apply custom flat
    with FlatModel(customflat) as flat:
        im,applied_flat = do_correction(im, flat)

    # add history entry
    time = datetime.now()
    stepdescription = 'Applied custom flat for NRCA5 detector; applyflat.py %s'%time.strftime('%Y-%m-%d %H:%M:%S')
    # writing to file doesn't save the time stamp or software dictionary 
    # with the History object, but left here for completeness
    software_dict = {'name':'applyflat.py',
                     'author':'Micaela Bagley',
                     'version':'1.0',
                     'homepage':'ceers.github.io'}
    substr = util.create_history_entry(stepdescription,
                                          software=software_dict)
    new_history = []
    for item in im.history:
        new_history.append(item)
    new_history.append(substr)
    im.history = new_history

    ## outputs
    # copy original
    log.info('Copying input to %s'%image.replace('_cal.fits', suffix))
    shutil.copy2(image, image.replace('_cal.fits', suffix))
    log.info('Saving output to %s'%image)
    im.save(image)
    fix_header(image)


def main():
    parser = argparse.ArgumentParser(description=
        'Apply a custom flat to NRCA5 cal images.')
    parser.add_argument('--output_dir', type=str, default='calibrated',
        help='Output directory. Default is calibrated.')
    parser.add_argument('--runone', type=str,
        help='Filename of single NRCA5 input file. If set, overrides the runall argument')
    parser.add_argument('--runall', action='store_true',
        help='Set to run all *nrca5_cal.fits images in the output_dir directory.')
    args = parser.parse_args()

    if args.runone:
        calimage = os.path.join(args.output_dir, args.runone)
        apply_custom_flat(calimage)

    elif args.runall:
        cals = glob(os.path.join(args.output_dir, '*_nrca5_cal.fits'))
        cals.sort()
        for calimage in cals:
            apply_custom_flat(calimage)
    
    else:
        print('Specify either --runone with a single input filename or --runall to process all *_nrca5_cal.fits files')


if __name__ == '__main__':
    main()


