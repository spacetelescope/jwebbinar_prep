##########################################################################
### JWST Pipeline practice - Reducing a subset of the CEERS 5 pointing ###
##########################################################################


We have created a Jupyter notebook (ceers_nircam_reduction.ipynb) that 
details how to reduce simulated CEERS NIRCam images through all stages of 
the JWST Calibration Pipeline. 

Jupyter notebook can be installed in your conda environment with:
    conda install -c conda-forge notebook
or 
    pip install notebook

You can then launch with:
    jupyter notebook

And click on ceers_nircam_reduction.ipynb from the file list. 
In the notebook, you will reduce 6 raw images (3 each of F115W and 
F277W imaging) and create 2 partial mosaics. 

***Note that running the steps in the notebook will take up ~15 GB of space!
   This total includes:
    - ~900 MB of raw images and auxiliary files (galaxy seed images)
    - ~400 MB of custom reference files
    - ~7-9 GB downloaded reference files in your CRDS_PATH directory
    - ~6.8 GB of pipeline outputs


############
### Set-up
See ceers_nircam_reduction.ipynb for more info.

1. Install the jwst calibration package, we recommend using version 1.3.3
   (replace <env_name> with your chosen environment name):

    conda create -n <env_name> python
    conda activate <env_name>
    pip install jwst==1.3.3

2. Set CRDS-related environment variables. For example, in bash:

    $ export CRDS_PATH=<your-path-to-crds-cache>/crds_cache
    $ export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
    $ export CRDS_CONTEXT='jwst_0764.pmap'


Raw images are in the uncals directory.
All pipeline outputs will be directed to the calibrated directory.


#####################
### Reduction steps
All steps, examples and explanations are available in the Jupyter notebook
ceers_nircam_reduction.ipynb


########################
### Directory contents

applyflat.py  - custom script to remove detector A5 feature, run between 
                Stage 2 (runimage2) and Stage 3 (runimage3)
calibrated  - directory for pipeline outputs
CEERSlogo.png  - CEERS logo image, embedded in Jupyter notebook
CEERSmap.png  - Image of CEERS observing configuration for June 2022, 
                embedded in Jupyter notebook
ceers_nircam_reduction.ipynb  - Jupyter notebook demonstrating reduction
                                steps
customflats  - directory of custom flat fields that are applied to A5 
               images to remove the feature present in simulated data
detector1_edited.asdf  - parameter file for Stage 1  (rundetector1)
f115w_nrca1.json  - association files grouping calibrated F115W images,
                    passed as inputs to Stage 3 to create partial F115W mosaic
f277w_nrca5.json - association files grouping calibrated F277W images,
                    passed as inputs to Stage 3 to create partial F277W mosaic
gains_v2.1.0  - directory of custom gain maps to replace the gain reference 
                files in CRDS
image2_edited.asdf  - parameter file for Stage 2 (runimage2)
image3_swc_edited.asdf  - parameter file for Stage 3 (runimage3) short 
                          wavelength filters 
image3_lwc_edited.asdf  - parameter file for Stage 3 (runimage3) long
                          wavelength filters 
jw0134500500*_01101_0000*_nrca*.json  - association files for running 
                                        SkyMatchStep on each calibrated image
                                        individually before Stage 3
plotimages.py  - convenience function to plot two images side-by-side for 
                 comparison of results before and after running a 
                 pipeline step
remstriping.py  - custom script to remove striping from countrate images,
                  run between Stage 1 and Stage 2 
rundithers  - batch script to run Stages 1 and 2 and custom steps on the 4
              remaining raw images that will be combined in the partial mosaics
skymatch_edited.asdf  - parameter file for SkyMatchStep
uncals  - directory of raw images and seed images
    *uncal.fits  - raw NIRCam images
    *seed_image.fits  - seed image from Mirage simulations, used to mask
                        source flux when measuring striping patterns


