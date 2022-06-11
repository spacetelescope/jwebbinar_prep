# JWebbinar material

Welcome to the maintained version of the JWebbinar notebooks.
If you are looking for the version that matches the recordings and PDF versions of the notebooks, you should navigate to the appropriate branch (branch "webbinar1" for JWebbinar 1, branch "webbinar2" for JWebbinar 2, and so on).

The majority of the notebooks in this repository run in an environment with the latest released versions of the JWST Calibration pipeline and Jdaviz (which installs astropy, specutils, matplotlib, and other packages required for JWST analysis).

To download the material, you can clone this repository:  
`git clone https://github.com/spacetelescope/jwebbinar_prep.git ~/my_jwebbinar_prep`  
`cd ~/my_jwebbinar_prep`  

Alternatively, you can manually download the individual notebooks by clicking on the pipeline_product_session folder, then click on each file, and right-click on "Raw" to select "Save link as".

To run these notebooks locally, you will need a Python environment set up that is capable of running Jupyter notebooks. If you are new to Python, you may want to look here: https://www.python.org/about/gettingstarted/ for some resources on how to install and learn the basics of Python. Depending on your preferences and system choices, you may find the install instructions there sufficient, but note that many scientists find it easier to use the Anaconda python distribution and package manager: https://www.anaconda.com/products/individual.

Using anaconda o miniconda, you first create an environment with the name of your choice (we use `jwebbinar` for this example):  
`conda create -n jwebbinar python=3.9`  
`conda activate jwebbinar`  
`pip install jwst`  
`pip install jdaviz`  

You can now run `jupyter notebook` in the folder where you have downloaded the webinar notebooks.

If you have problems or questions, feel free to reach out to the JWST Help Desk at http://jwsthelp.stsci.edu/.
